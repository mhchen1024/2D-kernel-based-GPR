function [delta_y,delta_y0,delta_u] = Reactor_Data_train(BatchNum, T, Ts, noiseSignalRatio)

    global Tj F 

    
    ntime=T/Ts;
    
    nrun=BatchNum;

    devrun=0.2;
    drun= randn(BatchNum,ntime)*devrun;
    tu=GBN(nrun,ntime,6,1,12,3);     
    y=zeros(nrun,ntime);
    u=zeros(nrun,ntime);
    load Reactor_Trajectory.mat;
    

    for i=1:1:nrun
        x0=[200 ; 0 ; 90];
        for j=1:1:ntime
            F = 5;  
            Tj=U(j)+tu(i,j);
            u(i,j)=Tj;
            [t,x] = ode45(@Reactor_model,[0 Ts],x0);
            x(size(x,1),:)';
            x0=x(size(x,1),:)';
            y(i,j)=x0(3,1);
        end
    end
    


%  plot(y(3,:))

    delta_y0 = y - ones(nrun,1)*X_nominal(:,3)';  
    delta_u = tu;
    
    stdV = std(reshape(drun, BatchNum*ntime,1));
    stdY0 = std(reshape(delta_y0, BatchNum*ntime,1));
%     delta_y = delta_y0 + drun/stdV*stdY0*sqrt(noiseSignalRatio);
    delta_y = delta_y0 + drun;
    
%     save('Reactor_data','delta_y','delta_y0','delta_u','nrun','ntime');
    
%     Y = delta_y;
%     U = delta_u;
%     t = 35;
%     Z_t   = Y(:,t);  
%     Phi_t = U(:,t:-1:1);  
%     g_t = inv(Phi_t'*Phi_t )*Phi_t'*Z_t;
%     var(Z_t - Phi_t*g_t)
%   

%     for t = 1:35
%         Z_t   = delta_y0(:,t);  
%         Phi_t = delta_u(:,t:-1:1);  
%         g_t = inv(Phi_t'*Phi_t )*Phi_t'*Z_t;
% %         G(1:t,t) = g_t;
%         G(t,t:-1:1) = g_t;
%         plot(g_t)
%         pause(0.2)
%         var(Z_t - Phi_t*g_t);
%     end


    
end