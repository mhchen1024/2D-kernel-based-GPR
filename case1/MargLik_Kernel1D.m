function f = MargLik_Kernel1D(alpha,Z_t,Phi_t,sigma,t)
% alpha(1) is c, alpha(2) is lambda, alpha(3) is rho, alpha(4)=beta
% global sigma  BatchNum Z_t Phi_t t
BatchNum = size(Z_t,1);
c = alpha(1);
lambda = alpha(2);
rho = alpha(3);

for m = 1:t
    for n = 1:t 
        P(m,n) = c*lambda^((m+n)/2)*rho^(abs(m-n));
    end
end

Cov = Phi_t*P*Phi_t' + sigma(t)*eye(BatchNum);
u   = chol(Cov);
tu  = inv(u);
B   = tu*tu';
Det = (diag(u).^2);

f = Z_t'*B*Z_t + sum(log(Det));
