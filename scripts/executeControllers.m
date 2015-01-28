
close all
clc

global Xk Uk K t

% X0 = funnel{7,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{8,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{14,1,2}.x(1,:)';% + [0.01;0.01;0];
X0 = funnel{1,1,2}.x(1,:)' + [0.1;0.1;0];

tt = [];  Xt = [];

kiter = 2;

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


%% 1

disp('computing trajectory...')
itrans = 1;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);

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
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);

for k = 1:length(t1)
    [isect1] = checkIntersection(tmpV1,tmpV,[],Xk1(k,:),eye(2));
    if isect1
        X0 = Xk1(k,:)';
        Xt = [Xt; Xk1(1:k,:)];
        tt = [tt; t1(1:k)];
        break
    end
end

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
count = 0;
while isempty(indx)
    kp = k + count;
    count = count + 1;
    kp
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk1(kp,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    
    % assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
    XinR
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
        end
    end
end
funIndx = indx(1);

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
fi2 = funIndx

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

if funIndx > 60
    disp('computing trajectory...')
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
        end
    end
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
    fi2(2) = funIndx
        
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
end

%% 3
disp('computing trajectory...')
itrans = 2;
[~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
    = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);

X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(2);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
fi3 = funIndx
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

kend = 14;  % sensor turns true

if funIndx > 40
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    n = itrans-1;
    indx = [];
    for m = 1:size(XbndR,1) % for all transition funnels
        if ~isempty(XbndR{m,n})
            indx = [indx;m];
        end
    end
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
    fi3(2) = funIndx
%     t03{2} = t;
%     Xk03{2} = Xk;
        
    X0 = Xk1(kend,:)';
    Xt = [Xt; Xk1(1:kend,:)];
    tt = [tt; t1(1:kend)+tt(end)];
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk1(kend,:)';
    Xt = [Xt; Xk1(1:kend,:)];
    tt = [tt; t1(1:kend)+tt(end)];
end

Xturn1 = Xk1(kend,:);

%% 4

disp('computing trajectory...')
X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

iReg = 2;
clear funIndx
% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
indx = [];
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
if ~isempty(indx)  % use inward funnel; else use transition
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
    fi4 = funIndx
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
    
    if funIndx > 60
        X0 = Xk1(end,:)';
        Xt = [Xt; Xk1];
        tt = [tt; t1+tt(end)];
        
        [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
        XbndR'
        for m = 1:length(XinR) % for all inward funnels
            if ~isempty(XinR{m})
                indx = [indx;m];
            end
        end
        funIndx = indx(1);
        [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
        fi4(2) = funIndx
        
        X0 = Xk1(end,:)';
        Xt = [Xt; Xk1(1:end,:)];
        tt = [tt; t1(1:end)+tt(end)];
        
        figure(1)
        plot(Xk1(:,1),Xk1(:,2))
        hold on
        plot(Xk(:,1),Xk(:,2),'r')
        drawnow
    else
        X0 = Xk1(end,:)';
        Xt = [Xt; Xk1(1:end,:)];
        tt = [tt; t1(1:end)+tt(end)];
    end
end


%% 5

disp('computing trajectory...')
itrans = 3;
[~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
    = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = itrans-1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
fi5 = funIndx
% t05 = t;
% Xk05 = Xk;

figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

kend = 90;  % sensor turns false
% kend = [];

if funIndx > 40
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    n = itrans-1;
    indx = [];
    for m = 1:size(XbndR,1) % for all transition funnels
        if ~isempty(XbndR{m,n})
            indx = [indx;m];
        end
    end
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
    fi5(2) = funIndx
%     t05{2} = t;
%     Xk05{2} = Xk;
        
    if isempty(kend)
        kend = length(t1);
    end
    X0 = Xk1(kend,:)';
    Xt = [Xt; Xk1(1:kend,:)];
    tt = [tt; t1(1:kend)+tt(end)];
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    if isempty(kend)
        kend = length(t1);
    end
    X0 = Xk1(kend,:)';
    Xt = [Xt; Xk1(1:kend,:)];
    tt = [tt; t1(1:kend)+tt(end)];
end

Xturn2 = Xk1(kend,:);

%% 6

disp('computing trajectory...')
X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

iReg = 2;
clear funIndx
% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
indx = [];
for m = 1:length(XinR) % for all inward funnels
    if ~isempty(XinR{m})
        indx = [indx;m];
    end
end
if ~isempty(indx)  % use inward funnel; else use transition
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
    fi6 = funIndx
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
    
    kend = 10;
    
    if funIndx > 60
        X0 = Xk1(kend,:)';
        Xt = [Xt; Xk1(1:kend,:)];
        tt = [tt; t1(1:kend)+tt(end)];
        
        [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
        XbndR'
        for m = 1:length(XinR) % for all inward funnels
            if ~isempty(XinR{m})
                indx = [indx;m];
            end
        end
        funIndx = indx(1);
        [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx);
        fi6(2) = funIndx
        
        X0 = Xk1(end,:)';
        Xt = [Xt; Xk1(1:end,:)];
        tt = [tt; t1(1:end)+tt(end)];
        
        figure(1)
        plot(Xk1(:,1),Xk1(:,2))
        hold on
        plot(Xk(:,1),Xk(:,2),'r')
        drawnow
    else
        X0 = Xk1(kend,:)';
        Xt = [Xt; Xk1(1:kend,:)];
        tt = [tt; t1(1:kend)+tt(end)];
    end
end



%% 7

disp('computing trajectory...')
itrans = 2;
[~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
    = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);

X0
Xk1(end,:)'
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk(kend,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = itrans-1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
fi7 = funIndx
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(1)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

if funIndx > 40
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
    
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
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx);
    fi7(2) = funIndx
%     t03{2} = t;
%     Xk03{2} = Xk;
        
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    figure(1)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
end


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


options.fill = 1;
options.color = [0.5 0.95 0.5];
options.shade = 0.5;
j = 1;
for i = 1:length(fi1);
    clear E
    for k = 1:length(funnel{fi1(i),j,2}.t)
        tmp = inv(funnel{fi1(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnel{fi1(i),j,2}.x(k,:)',tmp*funnel{fi1(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 1;
options.color = [0.9 0.3 0.3];
options.shade = 0.5;
j = 2;
for i = 1:length(fi2);
    clear E
    for k = 1:length(funnelIn{fi2(i),j,2}.t)
        tmp = inv(funnelIn{fi2(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnelIn{fi2(i),j,2}.x(k,:)',tmp*funnelIn{fi2(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 1;
options.color = [0.3 0.7 0.3];
options.shade = 0.5;
j = 2;
for i = 1:length(fi3);
    clear E
    for k = 1:length(funnel{fi3(i),j,2}.t)
        tmp = inv(funnel{fi3(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnel{fi3(i),j,2}.x(k,:)',tmp*funnel{fi3(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 1;
options.color = [0.7 0.3 0.3];
options.shade = 0.5;
j = 2;
for i = 1:length(fi4);
    clear E
    for k = 1:length(funnelIn{fi4(i),j,2}.t)
        tmp = inv(funnelIn{fi4(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnelIn{fi4(i),j,2}.x(k,:)',tmp*funnelIn{fi4(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 1;
options.color = [0.2 0.6 0.2];
options.shade = 0.5;
j = 2;
for i = 1:length(fi5);
    clear E
    for k = 1:length(funnel{fi5(i),j,2}.t)
        tmp = inv(funnel{fi5(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnel{fi5(i),j,2}.x(k,:)',tmp*funnel{fi5(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 0.5;
options.color = [0.3 0.7 0.3];
options.shade = 0.5;
j = 3;
for i = 1:length(fi5);
    clear E
    for k = 1:length(funnel{fi5(i),j,2}.t)
        tmp = inv(funnel{fi5(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnel{fi5(i),j,2}.x(k,:)',tmp*funnel{fi5(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 0.5;
options.color = [0.5 0.1 0.1];
options.shade = 0.5;
j = 2;
for i = 1:length(fi6);
    clear E
    for k = 1:length(funnelIn{fi6(i),j,2}.t)
        tmp = inv(funnelIn{fi6(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnelIn{fi6(i),j,2}.x(k,:)',tmp*funnelIn{fi6(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

options.fill = 0.5;
options.color = [0.5 0.95 0.5];
options.shade = 0.5;
j = 2;
for i = 1:length(fi7);
    clear E
    for k = 1:length(funnel{fi7(i),j,2}.t)
        tmp = inv(funnel{fi7(i),j,2}.P(:,:,k));
        tmp = (tmp+tmp')/2;
        E(k) = ellipsoid(funnel{fi7(i),j,2}.x(k,:)',tmp*funnel{fi7(i),j,2}.rho(k));
    end
    plotEllipse(projection(E,[1 0;0 1;0 0]),options.color);
end

% plot(Xk00(:,1),Xk00(:,2),'Color',[0.2 0.2 0.2],'LineWidth',2)
% plot(Xk01(:,1),Xk01(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk02(:,1),Xk02(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk03(:,1),Xk03(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk04(:,1),Xk04(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)
% plot(Xk05(:,1),Xk05(:,2),'Color',[0.6 0.6 0.6],'LineWidth',2)

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

plot(Xturn1(1),Xturn1(2),'y+','LineWidth',3,'MarkerSize',16)
if ~isempty(Xturn2), plot(Xturn2(1),Xturn2(2),'yx','LineWidth',3,'MarkerSize',16), end

set(gca, 'Box', 'off' )
% set(gca, 'TickDir', 'out')
