function f = MargLik_Kernel1DSS(alpha,Z_t,Phi_t,sigma,t)
% alpha(1) is c, alpha(2) is lambda
% global sigma N BatchNum Z_t Phi_t t

BatchNum = size(Z_t,1);
c = alpha(1);
lambda = alpha(2);

for m = 1:t
    for n = 1:t 
        P(m,n) = c*(0.5*lambda^(m+n+max(m,n)) - 1/6*lambda^(3*max(m,n)));
    end
end

Cov = Phi_t*P*Phi_t' + sigma(t)*eye(BatchNum);
u   = chol(Cov);
tu  = inv(u);
B   = tu*tu';
Det = (diag(u).^2);

f = Z_t'*B*Z_t + sum(log(Det));
