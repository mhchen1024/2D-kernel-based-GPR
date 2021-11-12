function f = MargLik_2D_DC(alpha,Big_Z,Psi,sigma,N,BatchNum)
% alpha(1) is c, alpha(2) is lambda, alpha(3) is rho, alpha(4)=beta
% global Big_Z Psi sigma N BatchNum

c = alpha(1);
lambda = alpha(2);
rho = alpha(3);
beta = alpha(4);


L = (1+N)*N/2;
for m = 1:L
    for n = 1:L 
        [i_m,j_m] = find_location(m,N);
        [i_n,j_n] = find_location(n,N);
        P(m,n) = c*lambda^((i_m+i_n)/2)*rho^(abs(i_m-i_n))*beta^(abs(j_m-j_n));
    end
end


Cov = Psi*P*Psi' + sigma*eye(BatchNum*N);
u   = chol(Cov);
tu  = inv(u);
B   = tu*tu';
% Det = (diag(u).^2);
% f = Big_Z'*B*Big_Z + sum(log(Det));
f = Big_Z'*B*Big_Z + 2*sum(log(diag(u)));



% %% QR
% u = chol(P,'lower');
% Rc = triu(qr([Rd(1:L+1,1:L)*u,Rd(1:L+1,L+1);sqrt(sigma)*eye(L),zeros(L,1)]));
% R1 = Rc(1:L,1:L);
% R2 = Rc(1:L,L+1);
% r = Rc(L+1,L+1);
% f = r^2/sigma+(N*BatchNum-L)*log(sigma)+2*real(log(det(R1)));
