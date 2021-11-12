function f = MargLik_2D_SS_M(alpha,Big_Z,Psi,N,BatchNum)
% alpha(1) is c, alpha(2) is lambda, alpha(3) is beta
% global Big_Z Psi sigma N BatchNum

c = alpha(1);
lambda = alpha(2);
beta = alpha(3);
sigma = alpha(4);
lengthscale = alpha(5);
% location = alpha(6);


L = (1+N)*N/2;
for m = 1:L
    for n = 1:L 
        [i_m,j_m] = find_location(m,N);
        [i_n,j_n] = find_location(n,N);
        P(m,n) = c*(0.5*lambda^(i_m+i_n+max(i_m,i_n)) - 1/6*lambda^(3*max(i_m,i_n)))*beta^(abs(j_m-j_n));
    end
end

for m = 1:(N*BatchNum)
    i_m = mod(m,BatchNum);
    for n = 1:(N*BatchNum)
        i_n = mod(n,BatchNum);
        if i_m == i_n
            j_m = floor((m-1)/BatchNum) + 1;
            j_n = floor((n-1)/BatchNum) + 1;
            K_e(m,n) = sigma*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
%             if j_m ~= j_n
%                 K_e(m,n) = sigma*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
%             else
%                 K_e(m,n) = sigma*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2)) + location*j_m;
%             end
% %             K_e(m,n) = sigma*(j_m-location)*(j_n-location)*exp(-(j_m-j_n)^2/(0.00001+2*lengthscale^2));
        end
    end
end

Cov = Psi*P*Psi' + K_e;
u   = chol(Cov);
tu  = inv(u);
B   = tu*tu';
Det = (diag(u).^2);

f = Big_Z'*B*Big_Z + sum(log(Det));
