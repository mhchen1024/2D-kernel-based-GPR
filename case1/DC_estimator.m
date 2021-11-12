function g_DC = DC_estimator(Y,U,sigma)
[BatchNum,N] = size(Y);
options = optimset('MaxFunEvals',2e5,'MaxIter',2e5);
g_DC = zeros(N,N); 
P_DC = zeros(N,N); 
x0 = [1,   0.80,0.00] ;
LB = [0,  0.72,-0.99];
UB = [1000,0.99,0.99];
for t = 1:N
    Z_t   = Y(:,t);
    Phi_t = U(:,t:-1:1);
    LikelihoodFun = @(eta)MargLik_DC(eta,Z_t,Phi_t,sigma,t);
    theta = fminsearchbnd(LikelihoodFun,x0,LB,UB,options);
    c = theta(1);
    lambda = theta(2);
    rho = theta(3);
    for i = 1:t
        for j = 1:t
            P_DC(i,j) = c*lambda^((i+j)/2)*rho^(abs(i-j));
        end
    end
    g_DC(1:t,t) = P_DC(1:t,1:t)*Phi_t'*inv(Phi_t*P_DC(1:t,1:t)*Phi_t' + sigma(t)*eye(BatchNum))*Z_t;
end