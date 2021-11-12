function f = MargLik_Inherit_DC(eta,Z_t,Phi_t,sigma,t,g_Inherit_DC,n)
% eta(1) is c   eta(2) is lambda,  eta(3) is rho
% global Z_t Phi_t sigma t BatchNum g_ap g_Inherit_DC n

BatchNum = size(Z_t,1);
c =  eta(1);
lambda = eta(2);
rho = eta(3);


for i = 1:t
    for j = 1:t
      P(i,j) = c*lambda^((i+j)/2)*rho^(abs(i-j));
    end
end

b0 = zeros(t,1);

if t < n+1
    Cov = (Phi_t*P*Phi_t' + sigma(t)*eye(BatchNum));
    u   = chol(Cov);
    tu  = inv(u);
    BB  = tu*tu';
    Det = (diag(u).^2);
    f = (Z_t)'*BB*(Z_t) + sum(log(Det)); 
else
    Theta_ap = zeros(t,n);
    for num = 1:n
        Theta_ap(1:t-num,num) = g_Inherit_DC(1:(t-num),t-num);
        Theta_ap(t-num+1:t,num) = g_Inherit_DC(t-num,t-num);
    end
    if n == 0
        b0 = zeros(t,1);
        S0 = diag(zeros(n,1));
    else
        b0 = Theta_ap*eta(1,3+1:3+n)';
        S0 = diag(eta(1,3+n+1:3+2*n)');
    end
    Cov = (Phi_t*(P+ Theta_ap*S0*Theta_ap')*Phi_t' + sigma(t)*eye(BatchNum));
    u   = chol(Cov);
    tu  = inv(u);
    BB  = tu*tu';
    Det = (diag(u).^2);
    f = (Z_t - Phi_t*b0)'*BB*(Z_t - Phi_t*b0) + sum(log(Det)); 
end
