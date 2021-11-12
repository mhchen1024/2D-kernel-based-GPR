function g_SS = SS_estimator(Y,U,sigma)
[BatchNum,N] = size(Y);
options = optimset('MaxFunEvals',2e5,'MaxIter',2e5);
g_SS = zeros(N,N); %脉冲响应矩阵
P_SS = zeros(N,N); %核函数
x0 = [1,  0.95] ;
LB = [0,  0.90];
UB = [1000,  1];
for t = 1:N
    Z_t   = Y(:,t);
    Phi_t = U(:,t:-1:1);
    LikelihoodFun = @(eta)MargLik_SS(eta,Z_t,Phi_t,sigma,t);
    % 超参数估计
    theta = fminsearchbnd(LikelihoodFun,x0,LB,UB,options);
    c = theta(1);
    lambda = theta(2);
    % 核矩阵
    for i = 1:t
        for j = 1:t
            P_SS(i,j) = c*(0.5*lambda^(i+j+max(i,j)) - 1/6*lambda^(3*max(i,j)));
        end
    end
    % 脉冲响应估计
    g_SS(1:t,t) = P_SS(1:t,1:t)*Phi_t'*inv(Phi_t*P_SS(1:t,1:t)*Phi_t' + sigma(t)*eye(BatchNum))*Z_t;
end