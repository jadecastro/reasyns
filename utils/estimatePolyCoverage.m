
kindx = 2;
imode = 4;


icover = [];
isect = [];
for itrans = 1:length(aut.trans)
    if aut.trans{itrans}(1) == imode
        itransMode = itrans;
        % for mm = 1:size(ellFunnel,1)
        % isect(mm,:) = isinternal_quick(ellFunnel{mm,itrans,1},qCover{itransReg}');
        % end
        icover = [icover; any(Isect{itrans,kindx},1)];
    end
end

Isect1{imode}(1:size(ellFunnelC,1),1:size(qCover{itransMode},1)) = 0;
for i = 1:size(ellFunnelC,1)
    i
    if all(~isempty(ellFunnelC{i,imode,kindx}))
        Isect1{imode}(i,:) = isinternal_quick(ellFunnelC{i,imode,kindx},qCover{itransMode}');
    end
end

insideIntersect = all(icover,1);
insideCentralRT = any(Isect1{imode},1);
Ncover1 = sum(~insideIntersect);  % Number of points not in the intersection region
icover1 = insideCentralRT | insideIntersect;  % Number of points covered by the inward-directed RT
%         icover1 = [all(icover,1); any(isect,1)];
NcoverAct = sum(icover1);

disp([' Percent of region covered: ',num2str(NcoverAct/Ncover),' , ',num2str(NcoverAct),' out of ',num2str(Ncover)]);
