% x3 = xMu(1,3);

global funnelI

clear Xp Yp Zp xx pp MM xx0 pp0 MM0 xxu ppu MMu xxg ppg MMg LHS LHS1
indx = 0;

% Evaluate the msspoly objects for plotting
% for k = 1:length(V)
%     x3 = funnelI.x(k,3);
%     XX = msubs(V',x(3),x3);
%     [xxtmp,pptmp,MMtmp] = decomp(XX(k));
%     xx{k} = xxtmp;
%     pp{k} = pptmp;
%     MM{k} = MMtmp;
% end
if ~isempty(X.mssInt{2})
    tmp = X.mssInt{2};
    if ~isempty(tmp)
        for k = 1:length(tmp)
            x3 = 0;
            XX(k) = msubs(tmp(k)',x(3),x3);
            [xxtmp,pptmp,MMtmp] = decomp(XX(k));
            xx1{k} = xxtmp;
            pp1{k} = pptmp;
            MM1{k} = MMtmp;
        end
    end
end

x1 = -1:0.1:15;
x2 = -1:0.1:15;
clear LHS LHS0 LHSu LHSb LHSg
for i = 1:length(x1)
    for j = 1:length(x2)
%         for k = 1:length(V)
%             LHS{k}(i,j) = MM{k}*(prod(repmat([x1(i) x2(j)],length(pp{k}),1).^pp{k},2));
%         end
        if ~isempty(X.mssInt{2})
            for k = 1:length(tmp)
                LHS1{k}(i,j) = MM1{k}*(prod(repmat([x1(i) x2(j)],length(pp1{k}),1).^pp1{k},2));
            end
        end
    end
end

% figure
axis equal
hold on
% for k = 1:length(V)
%     if rho_d(k) < 1e9
%         contour(x1,x2,LHS{k}'-rho_d(k),[0 0],'k')
%     end
% end
if ~isempty(X.mssInt{1})
    for k = 1:length(tmp)
        contour(x1,x2,LHS1{k}',[0 0],'r')
    end
end
