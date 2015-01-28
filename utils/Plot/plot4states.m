
% x3Vec = -pi:pi/3:pi;
xyzInit(3) = xyzInit(3) - 0;
x3val = sin(xyzInit(3));
x4val = cos(xyzInit(3));

for l1 = 1:Nloc
    if ~isempty(X0{l1})
        X0tmp = X0{l1};
        break
    end
end

clear Xp Yp Zp xx pp MM xx0 pp0 MM0 xxu ppu MMu xxg ppg MMg LHSbTmp
indx = 0;
if ~isempty(X{1})
    % Evaluate the msspoly objects for plotting
    XX = msubs(X{1}',x(3:4),[x3val;x4val]);
    for k = 1:length(XX)
        [xxtmp,pptmp,MMtmp] = decomp(XX(k));
        xx{k} = xxtmp;
        pp{k} = pptmp;
        MM{k} = MMtmp;
    end
end
X02 = msubs(X0tmp',x(3:4),[x3val;x4val]);
for k = 1:length(X02)
    [xxtmp,pptmp,MMtmp] = decomp(X02(k));
    xx0{k} = xxtmp;
    pp0{k} = pptmp;
    MM0{k} = MMtmp;
end
if ~isempty(Xue{1})
    XU = msubs(Xue{1}',x(3:4),[x3val;x4val]);
    for k = 1:length(Xue{1})
        [xxtmp,pptmp,MMtmp] = decomp(XU(k));
        xxu{k} = xxtmp;
        ppu{k} = pptmp;
        MMu{k} = MMtmp;
    end
end
XG = msubs(Xg{1}',x(3:4),[x3val;x4val]);
for k = 1:length(Xg{1})
    [xxtmp,pptmp,MMtmp] = decomp(XG(k));
    xxg{k} = xxtmp;
    ppg{k} = pptmp;
    MMg{k} = MMtmp;
end
for k = 1:length(barrier)
    barrier2 = msubs(barrier{k},x(3:4),[x3val;x4val]);
    [xxb{k},ppb{k},MMb{k}] = decomp(barrier2);
    [xxbTmp{k},ppbTmp{k},MMbTmp{k}] = decomp(barrier{k});
end

x1 = -10:0.5:10;
x2 = -10:0.5:10;
clear LHS LHS0 LHSu LHSb LHSg
for i = 1:length(x1)
    for j = 1:length(x2)
        if ~isempty(X{1})
            for k = 1:length(X{1})
                LHS{k}(i,j) = MM{k}*(prod(repmat([x1(i) x2(j)],length(pp{k}),1).^pp{k},2));
            end
        end
        for k = 1:length(X0tmp)
            LHS0{k}(i,j) = MM0{k}*(prod(repmat([x1(i) x2(j)],length(pp0{k}),1).^pp0{k},2));
        end
        if ~isempty(Xue{1})
            for k = 1:length(Xue{1})
                LHSu{k}(i,j) = MMu{k}*(prod(repmat([x1(i) x2(j)],length(ppu{k}),1).^ppu{k},2));
            end
        end
        for k = 1:length(barrier)
            LHSb{k}(i,j) = MMb{k}*(prod(repmat([x1(i) x2(j)],length(ppb{k}),1).^ppb{k},2));
        end
        Vx = xd - x1(i);
        Vy = yd - x2(j);
        % TODO: resolve each mode without hard-coding each of the limits
        
        if -(1/e*(-Vx*x3val + Vy*x4val) - wlim) < 0 && Nloc > 1
            LHSbTmp(j,i) = MMbTmp{2}*(prod(repmat([x1(i) x2(j) x3val x4val],length(ppbTmp{2}),1).^ppbTmp{2},2));
        elseif (1/e*(-Vx*x3val + Vy*x4val) + wlim) < 0 && Nloc > 2
            LHSbTmp(j,i) = MMbTmp{3}*(prod(repmat([x1(i) x2(j) x3val x4val],length(ppbTmp{3}),1).^ppbTmp{3},2));
        else
            LHSbTmp(j,i) = MMbTmp{1}*(prod(repmat([x1(i) x2(j) x3val x4val],length(ppbTmp{1}),1).^ppbTmp{1},2));
        end
        for k = 1:length(Xg{1})
            LHSg{k}(i,j) = MMg{k}*(prod(repmat([x1(i) x2(j)],length(ppg{k}),1).^ppg{k},2));
        end
    end
end

figure(10)
[~,H] = contourf(x1,x2,LHSbTmp,[0 0]);
set(get(H,'Children'),'FaceAlpha',0.25)

figure
axis equal
hold on
[~,H] = contourf(x1,x2,LHSbTmp,[0 0]);
set(get(H,'Children'),'FaceAlpha',0.25)
if ~isempty(X{1})
    for k = 1:length(X{1})
        contour(x1,x2,LHS{k}',[0 0],'k')
    end
end
for k = 1:length(X0tmp)
    contour(x1,x2,LHS0{k}',[0 0],'g','LineWidth',2)
end
if ~isempty(Xue{1})
    for k = 1:length(Xue{1})
        %     contour(x1,x2,lhsu{k}',0+eps,'r','linewidth',2)
        kp1 = mod(k,length(Xue{1}))+1;
        plot([reg.extVert(1,k) reg.extVert(1,kp1)],[reg.extVert(2,k) reg.extVert(2,kp1)],'r','linewidth',2)
    end
end
if ~isempty(Xui{1})
    for k = 1:length(Xui{1})
        %     contour(x1,x2,lhsu{k}',0+eps,'r','linewidth',2)
        kp1 = mod(k,length(Xui{1}))+1;
        plot([reg.intVert(1,k) reg.intVert(1,kp1)],[reg.intVert(2,k) reg.intVert(2,kp1)],'r','linewidth',2)
    end
end
for k = 1:length(barrier)
    contour(x1,x2,LHSb{k}',[0 0]+eps,'b:')
end
contour(x1,x2,LHSbTmp,[0 0],'m','LineWidth',2)

% for k = 1:length(Xg{1})
%     contour(x1,x2,LHSg{k}',0+eps,'r--')
% end
plot(xyzG(1),xyzG(2),'rx','LineWidth',2,'MarkerSize',10)
drawnow
LHSbTmp

%%

if plt3d
    % Evaluate the msspoly objects for plotting
    if ~isempty(X)
        for k = 1:length(X{1})
            [xxtmp,pptmp,MMtmp] = decomp(X{1}(k));
            xx{k} = xxtmp;
            pp{k} = pptmp;
            MM{k} = MMtmp;
        end
    end
    for k = 1:length(X0{1})
        [xxtmp,pptmp,MMtmp] = decomp(X0{1}(k));
        xx0{k} = xxtmp;
        pp0{k} = pptmp;
        MM0{k} = MMtmp;
    end
    for k = 1:length(Xue{1})
        [xxtmp,pptmp,MMtmp] = decomp(Xue{1}(k));
        xxu{k} = xxtmp;
        ppu{k} = pptmp;
        MMu{k} = MMtmp;
    end
    for k = 1:length(Xg{1})
        [xxtmp,pptmp,MMtmp] = decomp(Xg{1}(k));
        xxg{k} = xxtmp;
        ppg{k} = pptmp;
        MMg{k} = MMtmp;
    end
    
    clear xxb ppb MMb
    for i = 1:length(barrier)
        [xxb{i},ppb{i},MMb{i}] = decomp(barrier{i});
    end
    
    x1 = -4:0.1:4;
    x2 = -4:0.1:4;
    x3 = -pi:pi/5:pi;
    xs = sin(x3);
    xc = cos(x3);
    clear LHS LHS0 LHSu LHSb3 LHSg
    indx = 1;
    for i = 1:length(x1)
        for j = 1:length(x2)
            for m = 1:length(x3)
                %             for k = 1:length(X)
                %                 LHS{k}(i,j) = MM{k}*(prod(repmat([x1(i) x2(j)],length(pp{k}),1).^pp{k},2));
                %             end
                %             LHS0(i,j) = MM0*(prod(repmat([x1(i) x2(j)],length(pp0),1).^pp0,2));
                %             for k = 1:length(Xu)
                %                 LHSu{k}(i,j) = MMu{k}*(prod(repmat([x1(i) x2(j)],length(ppu{k}),1).^ppu{k},2));
                %             end
                
                if 1
                    Vx = xd - x1(i);
                    Vy = yd - x2(j);
                    if -(1/e*(-Vx*xs(m) + Vy*xc(m)) - wlim) < 0 && Nloc > 1
                        LHSb3(j,i,m) = MMb{2}*(prod(repmat([x1(i) x2(j) xs(m) xc(m)],length(ppb{2}),1).^ppb{2},2));
                    elseif (1/e*(-Vx*xs(m) + Vy*xc(m)) + wlim) < 0 && Nloc > 1
                        LHSb3(j,i,m) = MMb{3}*(prod(repmat([x1(i) x2(j) xs(m) xc(m)],length(ppb{3}),1).^ppb{3},2));
                    else
                        LHSb3(j,i,m) = MMb{1}*(prod(repmat([x1(i) x2(j) xs(m) xc(m)],length(ppb{1}),1).^ppb{1},2));
                    end
                else
                    LHSb3(j,i,m) = MMb{1}*(prod(repmat([x1(i) x2(j) xs(m) xc(m)],length(ppb{1}),1).^ppb{1},2));
                end
                %             for k = 1:length(Xg{1})
                %                 LHSg3{k}(j,i,m) = MMg{k}*(prod(repmat([x1(i) x2(j) x3(m)],length(ppg{k}),1).^ppg{k},2));
                %             end
            end
        end
    end
    
    % figure(10)
    % axis equal
    % hold on
    % for k = 1:length(X)
    %     contour(x1,x2,LHS{k}',1+eps)
    % end
    % contour(x1,x2,LHS0',1+eps)
    % for k = 1:length(Xu)
    %     contour(x1,x2,LHSu{k}',1+eps)
    % end
    
    figure
    isosurface(x1,x2,x3,LHSb3,0+eps)
    xlabel('x'), ylabel('y'), zlabel('\theta')
    %
    % figure
    % isosurface(x1,x2,x3,LHSg3{1},0+eps)
    % xlabel('x'), ylabel('y'), zlabel('\theta')
    % figure
    % isosurface(x1,x2,x3,LHSg3{2},0+eps)
    % xlabel('x'), ylabel('y'), zlabel('\theta')
    
end