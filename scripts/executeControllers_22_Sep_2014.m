
close all
clc

global Xk Uk K t

load inspectionTask_9_22_withoutshortcut_step2_fixed.mat

% aut.trans(19)=[];
% aut.trans(20)=[];

X0 = funnel{1,1}.x(1,:)' + [0.1;-0.2;0];

tt = [];  Xt = [];

figure(101)
clf
hold on
axis equal
set(gcf,'RendererMode','manual','Renderer','OpenGL');

vBndMin = min(vBnd{:});  vBndMax = max(vBnd{:});

set(gca,'FontSize',12)
xlabel('x_r'); ylabel('y_r'); zlabel('\theta')

color = {'g','b','y','b','y','b','y','b','y',[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1]};
Options.shade = 0.5;
hold on
for j = 1:length(pReg)
    Options.color = color{j};
    H=plot(pReg{j},Options);
    if any(j==1:9)
        set(H,'ZData',-10*ones(size(get(H,'XData'))))
    else
        set(H,'FaceAlpha',1)
    end
end

%% 1

itransList = [1:15 2:5];

kiter = 1;
disp('computing trajectory...')
iter = 0;
tend = 0;
foundReactiveX0 = [];
for itrans = itransList;
    iter = iter+1;
    % [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
    %     = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);
    [~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
        = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);
    
    for tidx = 19:length(aut.trans)
        aut.trans(length(aut.trans))=[];
    end
    X0
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    n = itrans;
%     indx = [];
%     for m = 1:size(XbndR,1) % for all transition funnels
%         if ~isempty(XbndR{m,n})
%             indx = [indx;m];
%         end
%     end
%     funIndx = indx(1);
    funIndx = 1;
    [t1{iter},Xk1{iter},t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
    
    X0 = Xk1{iter}(end,:)';
    Xt = [Xt; Xk1{iter}(1:end,:)];
    tt = [tt; t1{iter}(1:end)+tend];
    tend = tt(end);
    
    t01 = t;
    Xk01 = Xk;
    fi1 = funIndx
    
    if iter == 18
        % find a reactive funnel
        ireact = 16;
        [~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
            = getRegionsEllipses(kiter,ireact,aut,reg,funnel,[],funnelIn,[]);
        for tidx = 19:length(aut.trans)
            aut.trans(length(aut.trans))=[];
        end
        for trajidx = 1:length(t1{iter})
            [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(Xk1{iter}(trajidx,:)',regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
            for i = 1:length(XbndR)
                if ~isempty(XbndR{i})
                    break
                end
            end
        end
    end
    
    figure(1)
    plot(Xk1{iter}(:,1),Xk1{iter}(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
end

% foundReactiveX0 = Xk1{iter}(trajidx,:)';
foundReactiveX0 = funnel{1,16}.x(1,:)';

Xt(end-21:end,:) = [];
tt(end-21:end) = [];


%% 2 - reactive part

tt1 = [];  Xt1 = [];

itransList = [16:18];
X0 = foundReactiveX0;

kiter = 1;
disp('computing trajectory...')
iter = 0;
for itrans = itransList;
    iter = iter+1;
    % [~,~,~,~,~,~,~,reg1,reg2,~,~,Xbnd1,Xbnd2,~,~,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] ...
    %     = getRegionsEllipses(kiter,itrans,aut.trans,reg,funnel,ellFunnel,funnelIn,ellFunnelC);
    [~,~,~,~,~,~,~,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2] ...
        = getRegionsEllipses(kiter,itrans,aut,reg,funnel,[],funnelIn,[]);
    
    for tidx = 19:length(aut.trans)
        aut.trans(length(aut.trans))=[];
    end
    X0
    [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2);
    XbndR'
    n = itrans;
%     indx = [];
%     for m = 1:size(XbndR,1) % for all transition funnels
%         if ~isempty(XbndR{m,n})
%             indx = [indx;m];
%         end
%     end
%     funIndx = indx(1);
    funIndx = 1;
    [t1{iter},Xk1{iter},t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType);
    
    X0 = Xk1{iter}(end,:)';
    Xt = [Xt; Xk1{iter}(1:end,:)];
    tt = [tt; t1{iter}(1:end)+tend];
    Xt1 = [Xt1; Xk1{iter}(1:end,:)];
    tt1 = [tt1; t1{iter}(1:end)+tend];
    tend = tt(end);
    
    t01 = t;
    Xk01 = Xk;
    fi1 = funIndx
    
    figure(1)
    plot(Xk1{iter}(:,1),Xk1{iter}(:,2))
    hold on
    plot(Xk(:,1),Xk(:,2),'r')
    drawnow
end


%%
figure(101)
plot(Xt(:,1),Xt(:,2),'Color',[0.3 0.3 0.3],'LineWidth',2.5)
plot(Xt1(:,1),Xt1(:,2),'Color',[0.9 0 0],'LineWidth',2.5)

% set(gca, 'Box', 'off' )
% set(gca, 'TickDir', 'out')

ttmp = downsampleUniformly(tt,50)
for i = 1:length(ttmp)
    indx = find(tt == ttmp(i),1,'first');
    drawUnicycle(Xt(indx,:))
%     drawCar(Xt(indx,:))
end
