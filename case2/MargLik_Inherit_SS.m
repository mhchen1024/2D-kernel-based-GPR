function f = MargLik_Inherit_SS(eta,Z_t,Phi_t,sigma,t,g_Inherit_SS,n)

BatchNum = size(Z_t,1);

c = eta(1);
lambda = eta(2);

for i = 1:t
    for j = 1:t
      P(i,j) = c*(0.5*lambda^(i+j+max(i,j)) - 1/6*lambda^(3*max(i,j)));
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
        Theta_ap(1:t-num,num) = g_Inherit_SS(1:(t-num),t-num);
        Theta_ap(t-num+1:t,num) = g_Inherit_SS(t-num,t-num);
    end
    b0 = Theta_ap*eta(1,2+1:2+n)';
    S0 = diag(eta(1,2+n+1:2+2*n)');
    Cov = (Phi_t*(P+ Theta_ap*S0*Theta_ap')*Phi_t' + sigma(t)*eye(BatchNum));
    u   = chol(Cov);
    tu  = inv(u);
    BB  = tu*tu';
    Det = (diag(u).^2);
    f = (Z_t - Phi_t*b0)'*BB*(Z_t - Phi_t*b0) + sum(log(Det)); 
end
