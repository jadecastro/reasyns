function [vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2,...
    ellBnd11,ellBnd21,Vbnd1,Vbnd2,ellC11,ellC21,Vin1,Vin2] ... %,ellBnd11,ellBnd21,ellC11,ellC21] ...
    = getRegionsEllipses(km,itrans,aut,reg,funnel,ellFunnel,funnelIn,ellFunnelC)

[vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2] = deal([]);

count = 0;
for j = 1:length(reg)
    if ~(aut.q{aut.trans{itrans}(1)} == j) && ~(aut.q{aut.trans{itrans}(2)} == j)
        count = count+1;
        vAvoid{count} = reg{j}.v;
        vAvoidB{count} = reg{j}.vB;
        regAvoid{count} = reg{j};
    end
end
count = 0;
for j = 1:length(reg)
    if ~(aut.q{aut.trans{itrans}(1)} == j)
        count = count+1;
        vAvoid1{count} = reg{j}.v;
        vAvoid1B{count} = reg{j}.vB;
        regAvoid1{count} = reg{j};
    end
end
count = 0;
for j = 1:length(reg)
    if ~(aut.q{aut.trans{itrans}(2)} == j)
        count = count+1;
        vAvoid2{count} = reg{j}.v;
        vAvoid2B{count} = reg{j}.vB;
        regAvoid2{count} = reg{j};
    end
end

% NB: reach tube constructs are m x n cell arrays of funnels (arrays of p
% ellipsoids), where:
% m: # of funnels in the transition
% n: # of outgoing transitions for the current region
% p: # of ellipsoids making up the funnel

% Get outgoing reach tubes from initial state
count = 0;
[ellBnd1,ellBnd2,Vbnd1,Vbnd2,ellC1,ellC2,Vin1,Vin2,ellBnd11,ellBnd21,ellC11,ellC21] = deal([]);
for indx = 1:length(aut.trans)
    if aut.trans{itrans}(1) == aut.trans{indx}(1),
        count = count+1;
        count1 = 0;
        for ii = 1:size(funnel,1)
            if ~isempty(funnel{ii,indx,km})
                count1 = count1+1;
                for j = 1:length(funnel{ii,indx,km}.t)
                    tmp = inv(funnel{ii,indx,km}.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    ellBnd11{count1,count}(j,1) = ellipsoid(funnel{ii,indx,km}.x(j,:)',tmp*funnel{ii,indx,km}.rho(j))';
                end
%                 ellBnd1{count1,count} = ellFunnel{ii,indx,km};
                Vbnd1{count1,count} = ones(size(funnel{ii,indx,km}.V)) - funnel{ii,indx,km}.V;
            end
        end
    end
end
% Get outgoing reach tubes from final state
count = 0;
for indx = 1:length(aut.trans)
    if aut.trans{itrans}(2) == aut.trans{indx}(1),
        count = count+1;
        count1 = 0;
        for ii = 1:size(funnel,1)
            if ~isempty(funnel{ii,indx,km})
                count1 = count1+1;
                for j = 1:length(funnel{ii,indx,km}.t)
                    tmp = inv(funnel{ii,indx,km}.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    ellBnd21{count1,count}(j,1) = ellipsoid(funnel{ii,indx,km}.x(j,:)',tmp*funnel{ii,indx,km}.rho(j));
                end
%                 ellBnd2{count1,count} = ellFunnel{ii,indx,km};
                Vbnd2{count1,count} = ones(size(funnel{ii,indx,km}.V)) - funnel{ii,indx,km}.V;
            end
        end
    end
end
% Get inward-directed funnels for the two regions from the
% previous iteration
count1 = 0;
count2 = 0;
for ii = 1:size(funnelIn,1)
    if ~isempty(funnelIn{ii,aut.trans{itrans}(1),km})
        count1 = count1+1;
        for j = 1:length(funnelIn{ii,aut.trans{itrans}(1),km}.t)
            tmp = inv(funnelIn{ii,aut.trans{itrans}(1),km}.P(:,:,j));
            tmp = (tmp+tmp')/2;
            ellC11{count1}(j,1) = ellipsoid(funnelIn{ii,aut.trans{itrans}(1),km}.x(j,:)',tmp*funnelIn{ii,aut.trans{itrans}(1),km}.rho(j));
        end
%         ellC1{count1} = ellFunnelC{ii,trans{itrans}(1),km};
        Vin1{count1} = ones(size(funnelIn{ii,aut.trans{itrans}(1),km}.V)) - funnelIn{ii,aut.trans{itrans}(1),km}.V;
    end
    if ~isempty(funnelIn{ii,aut.trans{itrans}(2),km})
        count2 = count2+1;
        for j = 1:length(funnelIn{ii,aut.trans{itrans}(2),km}.t)
            tmp = inv(funnelIn{ii,aut.trans{itrans}(2),km}.P(:,:,j));
            tmp = (tmp+tmp')/2;
            ellC21{count2}(j,1) = ellipsoid(funnelIn{ii,aut.trans{itrans}(2),km}.x(j,:)',tmp*funnelIn{ii,aut.trans{itrans}(2),km}.rho(j));
        end
%         ellC2{count2} = ellFunnelC{ii,trans{itrans}(2),km};
        Vin2{count2} = ones(size(funnelIn{ii,aut.trans{itrans}(2),km}.V)) - funnelIn{ii,aut.trans{itrans}(2),km}.V;
    end
end