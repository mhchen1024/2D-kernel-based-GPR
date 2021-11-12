function Reactor_Trajectory_train(Ts,T)

    global Tj F

    ntime=T/Ts;

    for i=1:1:30
        YR(i,1)=90+(120-90)*i/(30);
    end
    for i=31:1:35
        YR(i,1)=120;
    end

    kp=1;
    nrun=20;
    U=zeros(ntime,1)+120;
    DU=zeros(ntime,1);
    Y=zeros(ntime,1);
    
    X_nominal = zeros(ntime,3);
    for ite = 1:20
        for i=1:1:nrun
            x0=[200 ; 0 ; 90];
            for k=1:1:ntime
                F=5;
                Tj=U(k);
                [t,x] = ode45(@Reactor_model,[0 Ts],x0);
                x0=x(size(x,1),:)';
                X_nominal(k,:) = x0';
                Y(k)=x0(3,1);
            end
            EK=YR-Y;    
            DU=kp*EK;
            U=U+DU;
        end   
    end
    save('Reactor_Trajectory','U','X_nominal');
% figure
% subplot(211);plot(Y(1:ntime),'linewidth',1);ylabel('r')
% subplot(212);plot(U(1:ntime),'linewidth',1);ylabel('T_J');xlabel('Sample')
