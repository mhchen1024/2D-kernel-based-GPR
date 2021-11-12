function f = MargLik_2D_SS(alpha,Big_Z,Psi,sigma,N,BatchNum)
% alpha(1) is c, alpha(2) is lambda, alpha(3) is beta
% global Big_Z Psi sigma N BatchNum

c = alpha(1);
lambda = alpha(2);
beta = alpha(3);

L = (1+N)*N/2;
for m = 1:L
    for n = 1:L 
        [i_m,j_m] = find_location(m,N);
        [i_n,j_n] = find_location(n,N);
        P(m,n) = c*(0.5*lambda^(i_m+i_n+max(i_m,i_n)) - 1/6*lambda^(3*max(i_m,i_n)))*beta^(abs(j_m-j_n));
    end
end

Cov = Psi*P*Psi' + sigma*eye(BatchNum*N);
u   = chol(Cov);
tu  = inv(u);
B   = tu*tu';
Det = (diag(u).^2);

f = Big_Z'*B*Big_Z + sum(log(Det));
