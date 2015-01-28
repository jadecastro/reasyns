
kindx = 2;
imode = 1;


icover = [];
isect = [];
for itrans = 1:length(aut.trans)
    if aut.trans{itrans}(1) == imode
        itransMode = aut.trans{itrans}(1);
        %                 for mm = 1:size(ellFunnel,1)
        %                     isect(mm,:) = isinternal(ellFunnel{mm,itrans,1},qCover{itransReg}');
        %                 end
        icover = [icover; any(Isect{itrans,1},1)];
    end
end

for i = 1:size(ellFunnelC,1)
    Isect1{imode}(i,:) = isinternal(ellFunnelC{i,imode,kindx},qCover{itransMode}');
end

insideIntersect = all(icover,1);
insideCentralRT = any(Isect1{imode},1);
Ncover1 = sum(~insideIntersect);  % Number of points not in the intersection region
icover1 = insideCentralRT.*(~insideIntersect(1:length(Isect1{imode})));  % Number of points covered by the inward-directed RT
%         icover1 = [all(icover,1); any(isect,1)];
NcoverAct = sum(icover1);

disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
