

clear Xp Yp Zp xx pp MM xx0 pp0 MM0 xxu ppu MMu xxg ppg MMg
indx = 0;
if ~isempty(X{1})
    % Evaluate the msspoly objects for plotting
    XX = X{1}';
    for k = 1:length(XX)
        [xxtmp,pptmp,MMtmp] = decomp(XX(k));
        xx{k} = xxtmp;
        pp{k} = pptmp;
        MM{k} = MMtmp;
    end
end
X02 = X0{2}';
for k = 1:length(X02)
    [xxtmp,pptmp,MMtmp] = decomp(X02(k));
    xx0{k} = xxtmp;
    pp0{k} = pptmp;
    MM0{k} = MMtmp;
end
XU = Xue{1}';
for k = 1:length(Xue{1})
    [xxtmp,pptmp,MMtmp] = decomp(XU(k));
    xxu{k} = xxtmp;
    ppu{k} = pptmp;
    MMu{k} = MMtmp;
end
XG = Xg{1}';
for k = 1:length(Xg{1})
    [xxtmp,pptmp,MMtmp] = decomp(XG(k));
    xxg{k} = xxtmp;
    ppg{k} = pptmp;
    MMg{k} = MMtmp;
end
for k = 1:length(barrier)
    barrier2 = barrier{k};
    [xxb{k},ppb{k},MMb{k}] = decomp(barrier2);
    [xxbTmp{k},ppbTmp{k},MMbTmp{k}] = decomp(barrier{k});
end

x1 = -1:0.01:1;
x2 = -1:0.01:1;
clear LHS LHS0 LHSu LHSb LHSg
for i = 1:length(x1)
    for j = 1:length(x2)
        if ~isempty(X{1})
            for k = 1:length(X{1})
                LHS{k}(i,j) = MM{k}*(prod(repmat([x1(i) x2(j)],length(pp{k}),1).^pp{k},2));
            end
        end
        for k = 1:length(X0{2})
            LHS0{k}(i,j) = MM0{k}*(prod(repmat([x1(i) x2(j)],length(pp0{k}),1).^pp0{k},2));
        end
        for k = 1:length(Xue{1})
            LHSu{k}(i,j) = MMu{k}*(prod(repmat([x1(i) x2(j)],length(ppu{k}),1).^ppu{k},2));
        end
        for k = 1:length(barrier)
            LHSb{k}(i,j) = MMb{k}*(prod(repmat([x1(i) x2(j)],length(ppb{k}),1).^ppb{k},2));
        end
        
        if x2(j) > 0.2 && Nloc > 1
            LHSbTmp(j,i) = MMbTmp{1}*(prod(repmat([x1(i) x2(j)],length(ppbTmp{1}),1).^ppbTmp{1},2));
        else
            LHSbTmp(j,i) = MMbTmp{2}*(prod(repmat([x1(i) x2(j)],length(ppbTmp{2}),1).^ppbTmp{2},2));
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
for k = 1:length(X0{2})
    contour(x1,x2,LHS0{k}',[0 0],'g','LineWidth',2)
end
for k = 1:length(Xue{1})
    %     contour(x1,x2,LHSu{k}',0+eps,'r','LineWidth',2)
    kp1 = mod(k,length(Xue{1}))+1;
    plot([reg.extVert(1,k) reg.extVert(1,kp1)],[reg.extVert(2,k) reg.extVert(2,kp1)],'r','LineWidth',2)
end
for k = 1:length(barrier)
    contour(x1,x2,LHSb{k}',[0 0]+eps,'b:')
end
contour(x1,x2,LHSbTmp,[0 0],'m','LineWidth',2)

drawnow
