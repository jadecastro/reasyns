% x3Vec = -pi:pi/3:pi;
x3Vec = 0;

clear Xp Yp Zp xx pp MM xx0 pp0 MM0 xxu ppu MMu xxg ppg MMg
indx = 0;
for x3i = x3Vec
    
    if ~isempty(X{1})
        % Evaluate the msspoly objects for plotting
        XX = msubs(X{1}',x(3),x3i);
        for k = 1:length(XX)
            [xxtmp,pptmp,MMtmp] = decomp(XX(k));
            xx{k} = xxtmp;
            pp{k} = pptmp;
            MM{k} = MMtmp;
        end
    end
    X02 = msubs(X0{1}',x(3),x3i);
    for k = 1:length(X02)
        [xxtmp,pptmp,MMtmp] = decomp(X02(k));
        xx0{k} = xxtmp;
        pp0{k} = pptmp;
        MM0{k} = MMtmp;
    end
    XU = msubs(Xu{1}',x(3),x3i);
    for k = 1:length(Xu{1})
        [xxtmp,pptmp,MMtmp] = decomp(XU(k));
        xxu{k} = xxtmp;
        ppu{k} = pptmp;
        MMu{k} = MMtmp;
    end
    XG = msubs(Xg{1}',x(3),x3i);
    for k = 1:length(Xg{1})
        [xxtmp,pptmp,MMtmp] = decomp(XG(k));
        xxg{k} = xxtmp;
        ppg{k} = pptmp;
        MMg{k} = MMtmp;
    end
    for k = 1:length(barrier)
        barrier2 = msubs(barrier{k},x(3),x3i);
        [xxb{k},ppb{k},MMb{k}] = decomp(barrier2);
    end
    
    x1 = -5:0.1:5;
    x2 = -5:0.1:5;
    clear LHS LHS0 LHSu LHSb LHSg
    for i = 1:length(x1)
        for j = 1:length(x2)
            if ~isempty(X{1})
                for k = 1:length(X{1})
                    LHS{k}(i,j) = MM{k}*(prod(repmat([x1(i) x2(j)],length(pp{k}),1).^pp{k},2));
                end
            end
            for k = 1:length(X0{1})
                LHS0{k}(i,j) = MM0{k}*(prod(repmat([x1(i) x2(j)],length(pp0{k}),1).^pp0{k},2));
            end
            for k = 1:length(Xu{1})
                LHSu{k}(i,j) = MMu{k}*(prod(repmat([x1(i) x2(j)],length(ppu{k}),1).^ppu{k},2));
            end
            for k = 1:length(barrier)
                LHSb{k}(i,j) = MMb{k}*(prod(repmat([x1(i) x2(j)],length(ppb{k}),1).^ppb{k},2));
            end
            for k = 1:length(Xg{1})
                LHSg{k}(i,j) = MMg{k}*(prod(repmat([x1(i) x2(j)],length(ppg{k}),1).^ppg{k},2));
            end
        end
    end
    
    figure
    axis equal
    hold on
    if ~isempty(X{1})
        for k = 1:length(X{1})
            contour(x1,x2,LHS{k}',0+eps,'k')
        end
    end
    for k = 1:length(X0{1})
        contour(x1,x2,LHS0{k}',0+eps,'g')
    end
    for k = 1:length(Xu{1})
        contour(x1,x2,LHSu{k}',0+eps,'r')
    end
    for k = 1:length(barrier)
        contour(x1,x2,LHSb{k}',0+eps,'b')
    end
    for k = 1:length(Xg{1})
        contour(x1,x2,LHSg{k}',0+eps,'r')
    end
    
end


%%

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
for k = 1:length(Xu{1})
    [xxtmp,pptmp,MMtmp] = decomp(Xu{1}(k));
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

x1 = -5:0.1:5;
x2 = -5:0.1:5;
x3 = -pi:pi/5:pi;
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
    if -(1/e*(-Vx*sinapprox(x3(m),n) + Vy*cosapprox(x3(m),n)) - 1) < 0
        LHSb3(j,i,m) = MMb{2}*(prod(repmat([x1(i) x2(j) x3(m)],length(ppb{2}),1).^ppb{2},2));
    elseif (1/e*(-Vx*sinapprox(x3(m),n) + Vy*cosapprox(x3(m),n)) - 1) < 0
        LHSb3(j,i,m) = MMb{3}*(prod(repmat([x1(i) x2(j) x3(m)],length(ppb{3}),1).^ppb{3},2));
    else
        LHSb3(j,i,m) = MMb{1}*(prod(repmat([x1(i) x2(j) x3(m)],length(ppb{1}),1).^ppb{1},2));
    end
else
    LHSb3(j,i,m) = MMb{1}*(prod(repmat([x1(i) x2(j) x3(m)],length(ppb{1}),1).^ppb{1},2));
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
