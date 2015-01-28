function isectBnd = chkContainment(funnelI,regAvoid,regSafe,regX,Hout,n,sampSkipValid)
% Checks to see if the funnel is misbehaving with respect to the polygonal
% avoid regions

isectBnd = false;

if sampSkipValid > 0
    chkIdx = 1:sampSkipValid:length(funnelI.t);
else
    chkIdx = length(funnelI.t):sampSkipValid:1;
end

for j = chkIdx
    tmp = inv(funnelI.P(:,:,j));
    tmp = (tmp+tmp')/2;
    E1 = ellipsoid(funnelI.x(j,:)',tmp*funnelI.rho(j));
    E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
    figure(5)
    plot(E)
    if any(intersect(E,regX.hExtB))
        isectBnd = true;
        break
    end
    for iii = 1:length(regAvoid)
        tmp1(iii) = size(extreme(regAvoid{iii}.p),1) == size(regAvoid{iii}.v,1);
    end
    if all(tmp1)  % avoid regions are convex
        for iii = 1:length(regAvoid)
            if any(intersect(E,regAvoid{iii}.pB2))
                isectBnd = true;
                break
            end
        end
    else  % they're not convex; take the safe regions (these should be convex)
        for iii = 1:length(regSafe)
            if any(intersect(E,regSafe{iii}.hBN2))
                isectBnd = true;
                break
            end
        end
    end
end

