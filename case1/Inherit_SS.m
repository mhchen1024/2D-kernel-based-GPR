function g_Inherit_SS = Inherit_SS(Y,U,sigma,n)

[BatchNum,N] = size(Y);
options = optimset('MaxFunEvals',2e5,'MaxIter',2e5);
g_Inherit_SS = zeros(N,N); 
g_ap     = zeros(N,1); 
P_Inherit_SS = zeros(N,N); 
ini = [1,0.95 ,zeros(1,2*n)];
lb  = [0.00, 0.0, -ones(1,n), zeros(1,n)];
ub  = [1000,1, ones(1,2*n)];
for t = 1:N
    Z_t   = Y(:,t);
    Phi_t = U(:,t:-1:1);
    
    LikelihoodFun = @(eta) MargLik_Inherit_SS(eta,Z_t,Phi_t,sigma,t,g_Inherit_SS,n);
    if t > n    
        [hyperparameter,fval,exitflag] = fminsearchbnd(LikelihoodFun,ini,lb,ub,options);
        b0 = hyperparameter(1,2+1:2+n)';
        S0 = diag(hyperparameter(1,2+n+1:2+2*n)');
        Theta_ap = zeros(t,n);
        for num = 1:n
            Theta_ap(1:t-num,num) = g_Inherit_SS(1:(t-num),t-num);
            Theta_ap(t-num+1:t,num) = g_Inherit_SS(t-num,t-num);
        end
        
        g_ap(1:t) = Theta_ap*b0;
    else
        [hyperparameter,fval,exitflag] = fminsearchbnd(LikelihoodFun,ini(1:2),lb(1:2),ub(1:2),options);
    end
    
    c = hyperparameter(1);
    lambda = hyperparameter(2);
    for i = 1:t
        for j = 1:t
            P_Inherit_SS(i,j) = c*(0.5*lambda^(i+j+max(i,j)) - 1/6*lambda^(3*max(i,j)));
        end
    end    
    
    if t > n
        
        Cov = (Phi_t*(P_Inherit_SS(1:t,1:t)+ Theta_ap*S0*Theta_ap')*Phi_t' + sigma(t)*eye(BatchNum));
        u   = chol(Cov);
        tu  = inv(u);
        BB  = tu*tu';
        Det = (diag(u).^2);
        
        beta = b0 + S0*Theta_ap'*Phi_t'*BB*(Z_t - Phi_t*g_ap(1:t));
        theta_t_r = P_Inherit_SS(1:t,1:t)*Phi_t'*BB*(Z_t - Phi_t*g_ap(1:t));
        g_Inherit_SS(1:t,t) = theta_t_r + Theta_ap*beta;
    else
        g_Inherit_SS(1:t,t) = P_Inherit_SS(1:t,1:t)*Phi_t'*inv(Phi_t*P_Inherit_SS(1:t,1:t)*Phi_t' + sigma(t)*eye(BatchNum))*(Z_t);      
    end
end