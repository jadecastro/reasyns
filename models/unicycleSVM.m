
% TODO: auto-generate SVM from symbolic math representation

global psi_dot k1 k2 k3 k4

Tfin = 100;

psi_dot = 0.001;

k1 = 1.1168;
k2 = -3.5294;
k3 = -1.5811;
k4 = -2.3677;


%% Linear model

A = [0 1 0 0;
    1/1.5 0 0 0;
    0 0 0 1;
    -13/3 0 0 0];
B = [0; -5/1.5; 0; 51/3];
C = eye(4);
D = zeros(4,1);
K = [k1 k3 k2 k4]; 

K = -lqr(A,B,1*eye(4),0.5);
k1 = K(1);  k3 = K(2);  k2 = K(3);  k4 = K(4);

%A = [0 1 0 0;0 0 1 0;0 0 0 1;-0.08 -0.06 -0.61 -0.2];
uniSys = ss(A,B,C,D);
uniCL = ss(A+B*K,zeros(4,1),C,D);

[yLin,tLin] = initial(uniCL,[0.4;0;0;0],Tfin);



%%
[t,y] = ode45(@uniModel,[0 Tfin],[0.4;0;0;0;0;0]);
%[t,y] = ode45(@uniModelLin,[0:0.001:Tfin],[0.004;0;0;0]);


figure(1)
subplot(211)
plot(tLin,yLin(:,1),'b',tLin,yLin(:,3),'r',t,y(:,1),'b--',t,y(:,5),'r--')
legend('\theta','\alpha')
subplot(212)
plot(tLin,yLin(:,2),'b',tLin,yLin(:,4),'r',t,y(:,2),'b--',t,y(:,6),'r--')
legend('\theta dot','\alpha dot')

figure(2)
subplot(211)
plot(t,y(:,1),t,y(:,3),t,y(:,5))
legend('\theta','\phi','\alpha')
subplot(212)
plot(t,y(:,2),t,y(:,4),t,y(:,6))
legend('\theta dot','\phi dot','\alpha dot')
