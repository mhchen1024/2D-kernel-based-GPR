function [Y,Y0,U,G] = generateData(batchNumber, sampleNumber, gbnTsw, noiseSignalRatio)


M = batchNumber;
N = sampleNumber;
ts = 1;
t_samples = (1:N)*ts;  % [1 2 3 4 5]
n = 1;


coefficient_K = [0.625];
coefficient_T = [0.001 -0.05 12];
val_K = polyval(coefficient_K, t_samples); 
val_T = polyval(coefficient_T, t_samples);
a1 = -exp(-ts./val_T);
b1 = val_K.*(1+a1);
trueLtvPara = [a1; b1]; 

gbnTsw = 5; 
U0 = gbngen(N*M, gbnTsw*ts);
U = reshape(U0, M, N);
Noise = randn(M,N);
Y0 = zeros(M,N);              

for k = 1:M  
    for t = 1:N  
        if t==1
            phi = [0, U(k,t)]; 
        else
            phi = [-Y0(k,t-1), U(k,t)];
        end
        theta = trueLtvPara(:,t);
        Y0(k,t) = phi*theta;  
    end
end

stdV = std(reshape(Noise, M*N,1));
stdY0 = std(reshape(Y0, M*N,1));
Y = Y0 + Noise/stdV*stdY0*sqrt(noiseSignalRatio); 


Y1 = zeros(M,N);           
for k = 1:M
    x = 0;
    for t = 1:N
        A(t) = [-a1(t)];
        B(t) = [b1(t)];
        C = [1];
        x = A(t)*x+B(t)*U(k,t);
        y = C*x;
        Y1(k,t) = y;
    end
end
G = zeros(N,N);
for t = 1:N
    for p = 1:t
        temp = 1;
        for q = 1:t-p
            temp = temp*A(t-q+1);
        end
        G(t,p) = C*temp*B(p);
    end   
end

