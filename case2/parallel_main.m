clc
clear all

core = 6;
p = parpool(core);

runs = 60;
F_inherit_DC = zeros(runs,1);
F_inherit_SS = zeros(runs,1);
F_LS = zeros(runs,1);
F_2D_DC_M = zeros(runs,1);
F_2D_SS_M = zeros(runs,1);
F_2D_DC_A = zeros(runs,1);
F_2D_SS_A = zeros(runs,1);

N = 35; 
BatchNum = 40;

options = optimset('MaxFunEvals',2e4,'MaxIter',2e4);
options = optimset('MaxFunEvals',2e4,'MaxIter',1e3);

G_inherit_DC = [];
G_inherit_SS = [];
G_LS = [];

G_2D_DC_M = [];
G_2D_DC_A = [];
G_2D_SS_M = [];
G_2D_SS_A = [];

error_list1 = [];
error_list2 = [];

noiseSignalRatio = 0.01;
gbnTsw=6;
T = 1.75;
Ts = 0.05
load truth_model.mat;
parfor iter = 1:runs
    disp(['第', num2str(iter),'轮仿真']);
    
    [Y, Y0, U] = Reactor_Data_train(BatchNum, T, Ts, noiseSignalRatio);
    

    Y = Y(:,1:N);
    Y0 = Y0(:,1:N);
    U = U(:,1:N);
    
    
%     G = zeros(N,N);
%     for t = 1:N
%         Z_t   = Y0(:,t);  % 列向量，包含所有批次t时刻的输出 
%         Phi_t = U(:,t:-1:1);  % 输入矩阵，每一行为 [u(t-1) u(t-2) ... u(0)]
%         g_t = inv(Phi_t'*Phi_t )*Phi_t'*Z_t;
% %         plot(g_t)
% %         pause(0.2)
%         G(t,t:-1:1) = g_t;
%     end
    
    g0 = [];
    for i = 1:N
        g0 = [g0;G(i,i:-1:1)'];
    end

    %% 1D kernel
    sigma = zeros(N,1);
    g_LS = zeros(N,N);
    for t = 1:N
        Z_t   = Y(:,t);
        Phi_t = U(:,t:-1:1);
        g_LS(1:t,t) = inv(Phi_t'*Phi_t)*Phi_t'*Z_t;

        sigma(t)    = var(Z_t - Phi_t*g_LS(1:t,t)); %sample variance
    end

    if sigma(N)>10^8 || sum(isnan(sigma))
        error_list1 = [error_list1;iter];
        continue;
    end
    
    %% DC kernel
   
    g_inherit_DC = [];
    g_inherit_SS = [];
    
    g_Inherit_DC_temp = Inherit_DC(Y,U,sigma,2);
    g_Inherit_SS_temp = Inherit_SS(Y,U,sigma,2);
    for i = 1:N
        g_inherit_DC = [g_inherit_DC;g_Inherit_DC_temp(1:i,i)];
        g_inherit_SS = [g_inherit_SS;g_Inherit_SS_temp(1:i,i)];
    end


    G_inherit_DC = [G_inherit_DC,g_inherit_DC];
    G_inherit_SS = [G_inherit_SS,g_inherit_SS];
    F_inherit_DC(iter) = 100*(1-(norm(g_inherit_DC-g0,2)/norm(g0,2)));
    F_inherit_SS(iter) = 100*(1-(norm(g_inherit_SS-g0,2)/norm(g0,2)));
    
    %% 2D 
    Big_Z = []; 
    Big_Z0 = [];  
    Psi = zeros(BatchNum*N,(N+1)*N/2);  
    for t = 1:N
        Big_Z = [Big_Z;Y(:,t)];
        Big_Z0 = [Big_Z0;Y0(:,t)];
        Psi((t-1)*BatchNum + 1: t*BatchNum,(t-1)*t/2+1:(t-1)*t/2+t) = U(:,t:-1:1);
    end
    theta_LS = inv(Psi'*Psi)*Psi'*Big_Z;  
    sigma = var(Big_Z - Psi*theta_LS); 
    if sigma>10^8 || sum(isnan(sigma))
        error_list2 = [error_list2;iter];
        continue;
    end
    

    
    %% 2D-DC-M kernel
    g_2D_DC_M = zeros((N+1)*N/2,1);
    P_2D_DC_M = zeros((N+1)*N/2,(N+1)*N/2);
    x0 = [1,0.80,0.00,0.00,0.5,0.00] ;
    LB = [0,0.72,-0.99,-0.99,0.00,0.00];
    UB = [1000,0.99,0.99,0.99,1,1];
    LikelihoodFun = @(eta) MargLik_2D_DC_M(eta,Big_Z,Psi,N,BatchNum);
    tic
    theta = fmincon(LikelihoodFun,x0,[],[],[],[],LB,UB);
    toc
    c = theta(1);
    lambda = theta(2);
    rho = theta(3);
    beta = theta(4);
    sigma_f = theta(5);
    lengthscale = theta(6);
    
    L = (1+N)*N/2;
    for m = 1:L
        for n = 1:L 
            [i_m,j_m] = find_location(m,N);
            [i_n,j_n] = find_location(n,N);
            P_2D_DC_M(m,n) = c*lambda^((i_m+i_n)/2)*rho^(abs(i_m-i_n))*beta^(abs(j_m-j_n));
        end
    end
    
    K_e = zeros(N*BatchNum,N*BatchNum);
    for m = 1:(N*BatchNum)
        i_m = mod(m,BatchNum);
        for n = 1:(N*BatchNum)
            i_n = mod(n,BatchNum);
            if i_m == i_n
                j_m = floor((m-1)/BatchNum) + 1;
                j_n = floor((n-1)/BatchNum) + 1;
                K_e(m,n) = sigma_f*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
            end
        end
    end
    
    g_2D_DC_M = P_2D_DC_M*Psi'*inv(Psi*P_2D_DC_M*Psi'+K_e)*Big_Z;
  
    %% 2D-DC-A kernel
    g_2D_DC_A = zeros((N+1)*N/2,1);
    P_2D_DC_A = zeros((N+1)*N/2,(N+1)*N/2);
    x0 = [1,0.72,0.99,0.00,0.5,0.00,1] ;
    LB = [0,0.72,-0.99,-0.99,0.00,0.00,0];
    UB = [1000,0.99,0.99,0.99,1,1,1000];
    LikelihoodFun = @(eta) MargLik_2D_DC_A(eta,Big_Z,Psi,N,BatchNum);
    tic
    theta = fmincon(LikelihoodFun,x0,[],[],[],[],LB,UB)
    toc
    c = theta(1);
    lambda = theta(2);
    rho = theta(3);
    beta = theta(4);
    sigma_f = theta(5);
    lengthscale = theta(6);
    c2 = theta(7);
    
    
    L = (1+N)*N/2;
    for m = 1:L
        for n = 1:L 
            [i_m,j_m] = find_location(m,N);
            [i_n,j_n] = find_location(n,N);
            P_2D_DC_A(m,n) = c*lambda^((i_m+i_n)/2)*rho^(abs(i_m-i_n)) + c2*beta^(abs(j_m-j_n));
        end
    end
    K_e = zeros(N*BatchNum,N*BatchNum);
    for m = 1:(N*BatchNum)
        i_m = mod(m,BatchNum);
        for n = 1:(N*BatchNum)
            i_n = mod(n,BatchNum);
            if i_m == i_n
                j_m = floor((m-1)/BatchNum) + 1;
                j_n = floor((n-1)/BatchNum) + 1;
                K_e(m,n) = sigma_f*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
            end
        end
    end
    
    g_2D_DC_A = P_2D_DC_A*Psi'*inv(Psi*P_2D_DC_A*Psi'+K_e)*Big_Z;


    
   
    %% 2D SS
    g_2D_SS_M = zeros((N+1)*N/2,1);
    P_2D_SS_M = zeros((N+1)*N/2,(N+1)*N/2);
    x0 = [1, 0.95,0.00,0.5,0.00];
    LB = [0, 0.90,-0.99,0.00,0.00];
    UB = [1000,1, 0.99,1,1];
    LikelihoodFun = @(eta) MargLik_2D_SS_M(eta,Big_Z,Psi,N,BatchNum);
    theta = fmincon(LikelihoodFun,x0,[],[],[],[],LB,UB)
    c = theta(1);
    lambda = theta(2);
    beta = theta(3);
    sigma_f = theta(4);
    lengthscale = theta(5);
    
    L = (1+N)*N/2;
    for m = 1:L
        for n = 1:L 
            [i_m,j_m] = find_location(m,N);
            [i_n,j_n] = find_location(n,N);
            P_2D_SS_M(m,n) = c*(0.5*lambda^(i_m+i_n+max(i_m,i_n)) - 1/6*lambda^(3*max(i_m,i_n)))*beta^(abs(j_m-j_n));
        end
    end
    K_e = zeros(N*BatchNum,N*BatchNum);
    for m = 1:(N*BatchNum)
        i_m = mod(m,BatchNum);
        for n = 1:(N*BatchNum)
            i_n = mod(n,BatchNum);
            if i_m == i_n
                j_m = floor((m-1)/BatchNum) + 1;
                j_n = floor((n-1)/BatchNum) + 1;
                K_e(m,n) = sigma*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
            end
        end
    end
    
    g_2D_SS_M = P_2D_SS_M*Psi'*inv(Psi*P_2D_SS_M*Psi'+K_e)*Big_Z;

    %% add
    g_2D_SS_A = zeros((N+1)*N/2,1);
    P_2D_SS_A = zeros((N+1)*N/2,(N+1)*N/2);
    x0 = [1, 0.95,0.98,0.5,0.00,0];
    LB = [0, 0.80,-0.99,0.00,0.00,0];
    UB = [1,1, 0.99,1,1,1];
    LikelihoodFun = @(eta) MargLik_2D_SS_A(eta,Big_Z,Psi,N,BatchNum);
    
    theta = fmincon(LikelihoodFun,x0,[],[],[],[],LB,UB)
    c = theta(1);
    lambda = theta(2);
    beta = theta(3);
    sigma_f = theta(4);
    lengthscale = theta(5);
    c2 = theta(6);
    
    L = (1+N)*N/2;
    for m = 1:L
        for n = 1:L 
            [i_m,j_m] = find_location(m,N);
            [i_n,j_n] = find_location(n,N);
            P_2D_SS_A(m,n) = c*(0.5*lambda^(i_m+i_n+max(i_m,i_n)) - 1/6*lambda^(3*max(i_m,i_n))) + c2*beta^(abs(j_m-j_n));
        end
    end
    K_e = zeros(N*BatchNum,N*BatchNum);
    for m = 1:(N*BatchNum)
        i_m = mod(m,BatchNum);
        for n = 1:(N*BatchNum)
            i_n = mod(n,BatchNum);
            if i_m == i_n
                j_m = floor((m-1)/BatchNum) + 1;
                j_n = floor((n-1)/BatchNum) + 1;
                K_e(m,n) = sigma_f*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
            end
        end
    end
    
    g_2D_SS_A = P_2D_SS_A*Psi'*inv(Psi*P_2D_SS_A*Psi'+K_e)*Big_Z;
    
    
    G_LS = [G_LS,theta_LS];
    G_2D_DC_M = [G_2D_DC_M,g_2D_DC_M];
    G_2D_DC_A = [G_2D_DC_A,g_2D_DC_A];
    G_2D_SS_M = [G_2D_SS_M,g_2D_SS_M];
    G_2D_SS_A = [G_2D_SS_A,g_2D_SS_A];
    F_LS(iter) = 100*(1-(norm(theta_LS-g0,2)/norm(g0,2)));
    F_2D_DC_M(iter) = 100*(1-(norm(g_2D_DC_M-g0,2)/norm(g0,2)));
    F_2D_SS_M(iter) = 100*(1-(norm(g_2D_SS_M-g0,2)/norm(g0,2)));
    F_2D_DC_A(iter) = 100*(1-(norm(g_2D_DC_A-g0,2)/norm(g0,2)));
    F_2D_SS_A(iter) = 100*(1-(norm(g_2D_SS_A-g0,2)/norm(g0,2)));

    
end
save('data_case2')
delete(p)