function u = GBN(nrun,ntime,tsw,tmin,tmax,ampl)

p = 1/tsw;
R = rand(nrun,ntime);

for i=1:1:nrun    
    if R(i,1) > 0.5 
        P_M = 1;
    else
        P_M = -1;
    end
    len = 1;
    for j = 1:1:ntime
        if R(i,j) < p 
            P_M = -P_M; 
        end
%    if ((R(i,j) < p) & (len >= tmin)) |(len >= tmax) 
%       P_M = -P_M; 
%       len = 1;
%    else
%       len = len+1;
%    end;
        u(i,j) = P_M*ampl;
    end
end