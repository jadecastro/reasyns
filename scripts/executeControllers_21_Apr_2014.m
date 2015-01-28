
close all
clc

global Xk Uk K t

load test_joined_23-Apr-2014.mat

funnelIn{1,1} = funnelJoin{3};
funnelIn{2,1} = funnelJoin{5};
funnelIn{3,2} = funnelJoin{1};
funnelIn{1,3} = funnelJoin{2};
funnelIn{1,4} = funnelJoin{4};

% X0 = funnel{7,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{8,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{14,1,2}.x(1,:)';% + [0.01;0.01;0];
X0 = funnelIn{1,1}.x(1,:)' + [0.1;-0.2;0];

tt = [];  Xt = [];

figure(1)
hold on

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'k--','LineWidth',2)
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'k--','LineWidth',2)
plot([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],'b:','LineWidth',2)
plot([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],'b:','LineWidth',2)
plot([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],'b:','LineWidth',2)
plot([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],'b:','LineWidth',2)

plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'k--','LineWidth',2)
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'k--','LineWidth',2)
plot([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],'b:','LineWidth',2)
plot([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],'b:','LineWidth',2)
plot([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],'b:','LineWidth',2)
plot([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],'b:','LineWidth',2)

set(gca,'FontSize',14)
xlabel('x_r'); ylabel('y_r')
text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,'r_1','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,'r_2','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,'r_3','FontSize',16)


kiter = 1;


%% 0

disp('computing trajectory...')
X0
iReg = 1;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = 1;%indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi0 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1];


%% 1

disp('computing trajectory...')
itrans = 1;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = itrans;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);

X0 = Xk1(end,:)';
Xt = [Xt; Xk1(1:end,:)];
tt = [tt; t1(1:end)+tt(end)];


t01 = t;
Xk01 = Xk;
fi1 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow
        
%% 2

disp('computing trajectory...')
X0
iReg = 2;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = 3;%indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi2 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 3
disp('computing trajectory...')
itrans = 2;
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];% clear nrm
% for i = 1:size(Xk,1), nrm(i) = norm(X0' - Xk(i,:)); end
% startIdx = find(min(nrm) == nrm)
% Xk(1:startIdx-1,:) = [];
% Uk(1:startIdx-1,:) = [];
% K(:,:,1:startIdx-1) = [];
% t(1:startIdx-1) = [];

    end
end
funIndx = 1;%indx(2);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
fi3 = funIndx;
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 4

disp('computing trajectory...')
X0
iReg = 3;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi4 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 5

disp('computing trajectory...')
itrans = 4;
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
fi5 = funIndx
% t05 = t;
% Xk05 = Xk;

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 6

disp('computing trajectory...')
X0
iReg = 4;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi6 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 7

disp('computing trajectory...')
itrans = 5;
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = 1;%indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
fi7 = funIndx
% t05 = t;
% Xk05 = Xk;

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];

%% 8

disp('computing trajectory...')
X0
iReg = 1;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = 2;

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi8 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];

%% 9

disp('computing trajectory...')
itrans = 1;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = itrans;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);

X0 = Xk1(end,:)';
Xt = [Xt; Xk1(1:end,:)];
tt = [tt; t1(1:end)];


t01 = t;
Xk01 = Xk;
fi9 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow
        
%% 10

disp('computing trajectory...')
X0
iReg = 2;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = 3;%indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi10 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 11
disp('computing trajectory...')
itrans = 2;
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);

X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];% clear nrm
% for i = 1:size(Xk,1), nrm(i) = norm(X0' - Xk(i,:)); end
% startIdx = find(min(nrm) == nrm)
% Xk(1:startIdx-1,:) = [];
% Uk(1:startIdx-1,:) = [];
% K(:,:,1:startIdx-1) = [];
% t(1:startIdx-1) = [];

    end
end
funIndx = 1;%indx(2);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
fi11 = funIndx;
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%% 12

disp('computing trajectory...')
X0
iReg = 3;
clear funIndx
indx = [];

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
funIndx = indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType);
fi12 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

X0 = Xk1(end,:)';
Xt = [Xt; Xk1];
tt = [tt; t1+tt(end)];


%%
figure
clf
hold on
axis equal

options.color = [0.9 0.9 0.9];
plot(pReg{1},options)
plot(pReg{3},options)

% options.fill = 1;
% options.color = [1 0.3 0.3];
% options.shade = 0.5;
% 
% j = 1;
% i = fi0;
% for k = 1:length(funnel{i,j,2}.t)
%     tmp = inv(funnel{i,j,2}.P(:,:,k));
%     tmp = (tmp+tmp')/2;
%     E(k) = ellipsoid(funnel{i,j,2}.x(k,:)',tmp*funnel{i,j,2}.rho(k));
% end
% plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);


%%
% options.fill = 1;
% options.color = [0.5 0.95 0.5];
% options.shade = 0.5;
% j = 1;
% for i = 1:length(fi1);
%     clear E
%     for k = 1:length(funnel{fi1(i),j,2}.t)
%         tmp = inv(funnel{fi1(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{fi1(i),j,2}.x(k,:)',tmp*funnel{fi1(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 1;
% options.color = [0.9 0.3 0.3];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi2);
%     clear E
%     for k = 1:length(funnelIn{fi2(i),j,2}.t)
%         tmp = inv(funnelIn{fi2(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnelIn{fi2(i),j,2}.x(k,:)',tmp*funnelIn{fi2(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 1;
% options.color = [0.3 0.7 0.3];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi3);
%     clear E
%     for k = 1:length(funnel{fi3(i),j,2}.t)
%         tmp = inv(funnel{fi3(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{fi3(i),j,2}.x(k,:)',tmp*funnel{fi3(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 1;
% options.color = [0.7 0.3 0.3];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi4);
%     clear E
%     for k = 1:length(funnelIn{fi4(i),j,2}.t)
%         tmp = inv(funnelIn{fi4(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnelIn{fi4(i),j,2}.x(k,:)',tmp*funnelIn{fi4(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 1;
% options.color = [0.2 0.6 0.2];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi5);
%     clear E
%     for k = 1:length(funnel{fi5(i),j,2}.t)
%         tmp = inv(funnel{fi5(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{fi5(i),j,2}.x(k,:)',tmp*funnel{fi5(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 0.5;
% options.color = [0.3 0.7 0.3];
% options.shade = 0.5;
% j = 3;
% for i = 1:length(fi5);
%     clear E
%     for k = 1:length(funnel{fi5(i),j,2}.t)
%         tmp = inv(funnel{fi5(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{fi5(i),j,2}.x(k,:)',tmp*funnel{fi5(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 0.5;
% options.color = [0.5 0.1 0.1];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi6);
%     clear E
%     for k = 1:length(funnelIn{fi6(i),j,2}.t)
%         tmp = inv(funnelIn{fi6(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnelIn{fi6(i),j,2}.x(k,:)',tmp*funnelIn{fi6(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end
% 
% options.fill = 0.5;
% options.color = [0.5 0.95 0.5];
% options.shade = 0.5;
% j = 2;
% for i = 1:length(fi7);
%     clear E
%     for k = 1:length(funnel{fi7(i),j,2}.t)
%         tmp = inv(funnel{fi7(i),j,2}.P(:,:,k));
%         tmp = (tmp+tmp')/2;
%         E(k) = ellipsoid(funnel{fi7(i),j,2}.x(k,:)',tmp*funnel{fi7(i),j,2}.rho(k));
%     end
%     plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
% end

% plot(Xk00(:,1),Xk00(:,2),'Color',[0.2 0.2 0.2],'LineWidth',2)
% plot(Xk01(:,1),Xk01(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk02(:,1),Xk02(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk03(:,1),Xk03(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk04(:,1),Xk04(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk05(:,1),Xk05(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)

%%
plot(Xt(:,1),Xt(:,2),'Color',[0.3 0.3 0.3],'LineWidth',2.5)


vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

plot([vReg{1}(2,1) vReg{1}(3,1)],[vReg{1}(2,2) vReg{1}(3,2)],'k--','LineWidth',2)
plot([vReg{1}(3,1) vReg{1}(4,1)],[vReg{1}(3,2) vReg{1}(4,2)],'k--','LineWidth',2)
plot([vRegB{1}(2,1) vRegB{1}(3,1)],[vRegB{1}(2,2) vRegB{1}(3,2)],'b:','LineWidth',2)
plot([vRegB{1}(3,1) vRegB{1}(4,1)],[vRegB{1}(3,2) vRegB{1}(4,2)],'b:','LineWidth',2)
plot([vRegBN{1}(2,1) vRegBN{1}(3,1)],[vRegBN{1}(2,2) vRegBN{1}(3,2)],'b:','LineWidth',2)
plot([vRegBN{1}(3,1) vRegBN{1}(4,1)],[vRegBN{1}(3,2) vRegBN{1}(4,2)],'b:','LineWidth',2)

plot([vReg{3}(1,1) vReg{3}(2,1)],[vReg{3}(1,2) vReg{3}(2,2)],'k--','LineWidth',2)
plot([vReg{3}(4,1) vReg{3}(1,1)],[vReg{3}(4,2) vReg{3}(1,2)],'k--','LineWidth',2)
plot([vRegB{3}(1,1) vRegB{3}(2,1)],[vRegB{3}(1,2) vRegB{3}(2,2)],'b:','LineWidth',2)
plot([vRegB{3}(4,1) vRegB{3}(1,1)],[vRegB{3}(4,2) vRegB{3}(1,2)],'b:','LineWidth',2)
plot([vRegBN{3}(1,1) vRegBN{3}(2,1)],[vRegBN{3}(1,2) vRegBN{3}(2,2)],'b:','LineWidth',2)
plot([vRegBN{3}(4,1) vRegBN{3}(1,1)],[vRegBN{3}(4,2) vRegBN{3}(1,2)],'b:','LineWidth',2)

set(gca,'FontSize',14)
xlabel('x_r'); ylabel('y_r')
text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,'r_1','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,'r_2','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,'r_3','FontSize',16)

plot(Xt(1,1),Xt(1,2),'ys','LineWidth',3,'MarkerSize',16)

set(gca, 'Box', 'off' )
% set(gca, 'TickDir', 'out')
