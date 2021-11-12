
function U = gbngen(N,Tsw,Seeds)

%%% gbngen.m; Matlab M file for generating GBN signal %%%
%
%    U = prbs(N,Tsw,Seeds)
%
%  This function generates GBN signal U in [-1 1] with Tmin = 1
%
%  Input arguments:
%   N :   Number of samples
%   Tsw :   Average switching time in samples 
%   Seeds : Seeds for random generator
%

if nargin < 2
   error('Not enough input arguments, try again')
end

% Detemine switching probability
psw = 1/Tsw;

if nargin > 2 
  rand('seed',Seeds); 
end;  

R = rand(N,1);

% Determine the initial value
if R(1) > 0.5, 
  P_M = 1;
else
  P_M = -1;
end;

U = zeros(N,1);

for k = 1:N,
   if (R(k) < psw)
      P_M = -P_M; 
   end;
   U(k) = P_M;
end

