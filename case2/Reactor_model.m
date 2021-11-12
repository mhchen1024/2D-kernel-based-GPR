% Fed-Batch  Reactor

function xdot=Reactor_model(t,x)

global Tj F

Vr= x(1,1);
z= x(2,1);
Tr= x(3,1);

Dr=2;
A=Vr*4/Dr;

% Parameters:
rho=50;
M=50;
E=30000;
cp=0.75;
lam=-20000;
U=100;
alpha=4.354e11;
R=1.986;
k=alpha*exp(-E/(R*(Tr+459.67)));

z0=1;
T0=90;

% Compute xdot:
xdot(1,1) = F;
xdot(2,1) = F*(z0-z)/Vr-k*z;
xdot(3,1) = F*(T0-Tr)/Vr-lam*k*z/(M*cp)-U*A*(Tr-Tj)/(rho*cp*Vr);
