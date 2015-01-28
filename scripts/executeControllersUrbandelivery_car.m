
close all
clc
clear all

load antsy4_urbanDelivery_car

global Xk Uk K t Xss Uss Kss

% X0 = funnel{7,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{8,1,2}.x(1,:)';% + [0.01;0.01;0];
% X0 = funnel{14,1,2}.x(1,:)';% + [0.01;0.01;0];
X0 = funnel{2,1,1}.x(1,:)' + [0.1;0.1;0;0];

tt = [];  Xt = [];

% kiter = 2;
kiter = 1;

plotWS(pReg)

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});
plot([vBndMin(1) vBndMin(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMin(2) vBndMin(2)],'k--','LineWidth',2)
plot([vBndMax(1) vBndMax(1)],[vBndMin(2) vBndMax(2)],'k--','LineWidth',2)
plot([vBndMin(1) vBndMax(1)],[vBndMax(2) vBndMax(2)],'k--','LineWidth',2)

for i = setdiff(1:length(vReg),4)
    for j = 1:size(vReg{i},1)
        plot([vReg{i}(j,1) vReg{i}(mod(j,size(vReg{i},1))+1,1)],[vReg{i}(j,2) vReg{i}(mod(j,size(vReg{i},1))+1,2)],'k--','LineWidth',2)
    end
end
set(gca,'FontSize',14)
xlabel('x_r'); ylabel('y_r')
text(vReg{1}(1,1)+0.2,vReg{1}(1,2)+0.5,'S','FontSize',16)
text(vReg{1}(2,1)+0.2,vReg{1}(1,2)+0.5,'O','FontSize',16)
text(vReg{2}(1,1)+0.2,vReg{2}(1,2)+0.5,'R','FontSize',16)
text(vReg{3}(1,1)+0.2,vReg{3}(1,2)+0.5,'C','FontSize',16)


%% 1

disp('computing trajectory...')
itrans = 1;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

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

for k = 1:length(t1)
    [isect1] = checkIntersection(tmpV1,tmpV,[],Xk1(k,:),eye(2));
    if isect1
        X01 = Xk1(k,:)';
        X0 = Xk(k,:)';
        Xt = [Xt; Xk1(1:k,:)];
        tt = [tt; t1(1:k)];
        break
    end
end

t01 = t;
Xk01 = Xk;
fi1 = funIndx

figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow
drawCar(Xk1(1,:))

        
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

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
fi2 = funIndx

figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

if funIndx > 1000 % optional depth 
    disp('computing trajectory...')
    X01 = Xk1(end,:)';
    X0 = Xk(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X01,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
        end
    end
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
    fi2(2) = funIndx
        
    X01 = Xk1(end,:)';
    X0 = Xk(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk(end,:)'; % NB: using last point in funnel trajectory
    X01 = Xk1(end,:)'; 
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
end

%% 3
disp('computing trajectory...')
itrans = 2;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

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
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
fi3 = funIndx
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

for k = 1:length(t1)
    [isect1] = checkIntersection(tmpV1,tmpV,[],Xk1(k,:),eye(2));
    if isect1
        X01 = Xk1(k,:)';
        X0 = Xk(k,:)';
        Xt = [Xt; Xk1(1:k,:)];
        tt = [tt; t1(1:k)+tt(end)];
        break
    end
end

kend = 14;  % sensor turns true

% if funIndx > 1000
%     X01 = Xk1(end,:)';
%     X0 = Xk(end,:)';
%     Xt = [Xt; Xk1];
%     tt = [tt; t1+tt(end)];
%     
%     [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
%     XbndR'
%     n = itrans-1;
%     indx = [];
%     for m = 1:size(XbndR,1) % for all transition funnels
%         if ~isempty(XbndR{m,n})
%             indx = [indx;m];
%         end
%     end
%     funIndx = indx(1);
%     [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
%     fi3(2) = funIndx
% %     t03{2} = t;
% %     Xk03{2} = Xk;
%         
%     X01 = Xk1(end,:)';
%     X0 = Xk(kend,:)';
%     Xt = [Xt; Xk1(1:kend,:)];
%     tt = [tt; t1(1:kend)+tt(end)];
%     
%     figure(4)
%     plot(Xk1(:,1),Xk1(:,2))
%     hold on
%     plot(Xk(:,1),Xk(:,2),'r')
%     drawnow
% else
%     X01 = Xk1(kend,:)';
%     X0 = Xk(kend,:)';
%     Xt = [Xt; Xk1(1:kend,:)];
%     tt = [tt; t1(1:kend)+tt(end)];
% end

% Xturn1 = Xk1(kend,:);

%% 4

disp('computing trajectory...')
X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X01,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

iReg = 3;
clear funIndx
% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
indx = [];
funType = [];
count = 0;
while isempty(indx)
    kp = k + count;
    count = count + 1;
    kp
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk1(kp,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    
    % assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
    XinR
    XbndR'
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
            funType = 1;
        end
    end
    if isempty(indx)
        for m = 1:length(XbndR) % for all transition funnels
            if ~isempty(XbndR{m})
                indx = [indx;m];
                funType = 2;
            end
        end
    end
end
funType
funIndx = indx(1);

if funType == 1
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
else
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnel,funIndx,modelType);
end
fi4 = funIndx

figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on
plot(Xk(:,1),Xk(:,2),'r')
drawnow

if funIndx > 1000
    X01 = Xk1(end,:)';
    X0 = Xk(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X01,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
        end
    end
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
    fi4(2) = funIndx
    
    X01 = Xk1(end,:)';
    X0 = Xk(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X01 = Xk1(end,:)';
    X0 = Xk(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
end

Xksav = Xk;
Xk1sav = Xk1;



%% 5

disp('computing trajectory...')
itrans = 3;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(2,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk(end,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
% funIndx = indx(1);
funIndx = 8; % FORCING THIS!!
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
fi5 = funIndx
% t05 = t;
% Xk05 = Xk;

figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

for k = 1:length(t1)
    [isect1] = checkIntersection(tmpV1,tmpV,[],Xk1(k,:),eye(2));
    if isect1
        X01 = Xk1(k,:)';
        X0 = Xk(k,:)';
        break
    end
end

kend = 90;  % sensor turns false
% kend = [];

% if funIndx > 1000
%     X01 = Xk1(end,:)';
%     X0 = Xk(end,:)';
%     Xt = [Xt; Xk1];
%     tt = [tt; t1+tt(end)];
%     
%     [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
%     XbndR'
%     n = itrans-1;
%     indx = [];
%     for m = 1:size(XbndR,1) % for all transition funnels
%         if ~isempty(XbndR{m,n})
%             indx = [indx;m];
%         end
%     end
%     funIndx = indx(1);
%     [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
%     fi5(2) = funIndx
% %     t05{2} = t;
% %     Xk05{2} = Xk;
%         
%     if isempty(kend)
%         kend = length(t1);
%     end
%     X01 = Xk1(kend,:)';
%     X0 = Xk(kend,:)';
%     Xt = [Xt; Xk1(1:kend,:)];
%     tt = [tt; t1(1:kend)+tt(end)];
%     
%     figure(4)
%     plot(Xk1(:,1),Xk1(:,2))
%     hold on
%     plot(Xk(:,1),Xk(:,2),'r')
%     drawnow
% else
%     if isempty(kend)
%         kend = length(t1);
%     end
%     X0 = Xk(kend,:)';
%     X01 = Xk1(kend,:)';
%     Xt = [Xt; Xk1(1:kend,:)];
%     tt = [tt; t1(1:kend)+tt(end)];
% end

% Xturn2 = Xk1(kend,:);

Xksav = Xk;
Xk1sav = Xk1;



%% 6
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

disp('computing trajectory...')
X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

iReg = 4;
clear funIndx
% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
indx = [];
count = 0;
% while isempty(indx)
%     kp = k + count;
%     count = count + 1;
%     kp
%     [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk(kp,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
%     % TODO: make sure we're in the right region!!
% 
%     % assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
%     XinR
%     for m = 1:length(XinR) % for all inward funnels
%         if ~isempty(XinR{m})
%             indx = [indx;m];
%         end
%     end
% end

% funIndx = indx(1);
funIndx = 18; % FORCING THIS!!
kp = 50;
Xt = [Xt; Xk1(1:kp,:)];
tt = [tt; t1(1:kp)+tt(end)];

[t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(Xk1(kp,:)',iReg,kiter,funnelIn,funIndx,modelType);
fi6 = funIndx

figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on
plot(Xk(:,1),Xk(:,2),'r')
drawnow

kend = 10;

if funIndx > 1000
    X0 = Xk(kend,:)';
    X01 = Xk1(kend,:)';
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
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
    fi6(2) = funIndx
    
    X0 = Xk(end,:)';
    X01 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    %         X0 = Xk(kend,:)';
    %         X01 = Xk1(kend,:)';
    %         Xt = [Xt; Xk1(1:kend,:)];
    %         tt = [tt; t1(1:kend)+tt(end)];
    X0 = Xk(end,:)';
    X01 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
end

Xksav1 = Xk;
Xk1sav1 = Xk1;

return

%% 7

disp('computing trajectory...')
itrans = 4;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

X0
Xk1(end,:)'
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk(end,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 1;  % go to state 1
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(1);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
fi7 = funIndx
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

kend = 30;

if funIndx > 1000
    X0 = Xk(end,:)';
    X01 = Xk1(end,:)';
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
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
    fi7(2) = funIndx
%     t03{2} = t;
%     Xk03{2} = Xk;
        
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk(kend,:)';
    X01 = Xk1(kend,:)';
    Xt = [Xt; Xk1(1:kend,:)];
    tt = [tt; t1(1:kend)+tt(end)];
    X01all = Xk1(end,:)';
    k = kend;  % force k in this funnel to be the time of the switch
end

Xturn1 = Xk1(kend,:)';
return 

%% 8

disp('computing trajectory...')
X0
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);

iReg = 4;
clear funIndx
% assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
XinR
indx = [];
count = 0;
while isempty(indx)
    kp = k + count;
    count = count + 1;
    kp
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk(kp,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    
    % assume we will not find a transition funnel as soon as we cross the border, so we look for inward funnels
    XinR
    for m = 1:length(XinR) % for all inward funnels
        if ~isempty(XinR{m})
            indx = [indx;m];
        end
    end
end

if ~isempty(indx)  % use inward funnel; else use transition
    funIndx = indx(1);
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
    fi8 = funIndx
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
    
    kend = 10;
    
    if funIndx > 1000
        X0 = Xk(kend,:)';
        X01 = Xk1(kend,:)';
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
        [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X01,iReg,kiter,funnelIn,funIndx,modelType);
        fi8(2) = funIndx
        
        X0 = Xk(end,:)';
        X01 = Xk1(end,:)';
        Xt = [Xt; Xk1(1:end,:)];
        tt = [tt; t1(1:end)+tt(end)];
        
        figure(4)
        plot(Xk1(:,1),Xk1(:,2))
        hold on
        plot(Xk(:,1),Xk(:,2),'r')
        drawnow
    else
%         X0 = Xk(kend,:)';
%         X01 = Xk1(kend,:)';
%         Xt = [Xt; Xk1(1:kend,:)];
%         tt = [tt; t1(1:kend)+tt(end)];
        X0 = Xk(end,:)';
        X01 = Xk1(end,:)';
        Xt = [Xt; Xk1(1:end,:)];
        tt = [tt; t1(1:end)+tt(end)];
    end
end

X0sav = Xk(end,:)';
X01sav = Xk1(end,:)';


%% 9

disp('computing trajectory...')
itrans = 5;
% [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
%     = getRegionsEllipses(2,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);
[~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
    = getRegionsEllipses(kiter,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelIn);

X0
Xk1(end,:)'
[XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
XbndR'
n = 2;  % go to state 5
indx = [];
for m = 1:size(XbndR,1) % for all transition funnels
    if ~isempty(XbndR{m,n})
        indx = [indx;m];
    end
end
funIndx = indx(2);
[t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
fi9 = funIndx
% t03{1} = t;
% Xk03{1} = Xk;
    
figure(4)
plot(Xk1(:,1),Xk1(:,2))
hold on 
plot(Xk(:,1),Xk(:,2),'r')
drawnow

kend = 20;

if funIndx > 1000
    X0 = Xk(end,:)';
    X01 = Xk1(end,:)';
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
    [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X01,itrans,kiter,funnel,funIndx,modelType);
    fi9(2) = funIndx
%     t03{2} = t;
%     Xk03{2} = Xk;
        
    X0 = Xk1(end,:)';
    Xt = [Xt; Xk1];
    tt = [tt; t1+tt(end)];
    
    figure(4)
    plot(Xk1(:,1),Xk1(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
else
    X0 = Xk(end,:)';
    X01 = Xk1(end,:)';
    Xt = [Xt; Xk1(1:end,:)];
    tt = [tt; t1(1:end)+tt(end)];
end

return 

%%

figure(5)

plot(Xt(:,1),Xt(:,2),'Color',[0.3 0.3 0.3],'LineWidth',2.5)
ttmp = downsampleUniformly(tt,30)
for i = 1:length(ttmp)
    indx = find(tt == ttmp(i),1,'first');
%     drawUnicycle(Xt(indx,:))
    drawCar(Xt(indx,:))
end
plot(Xturn1(1),Xturn1(2),'ro','LineWidth',3,'MarkerSize',10)

%%
figure(5)
clf
hold on
axis equal

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
