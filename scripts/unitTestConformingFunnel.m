

global ME

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = 5; % options.maxTrials2;
trans = vertcat(aut.trans{:});
options

IsectInAll{iModeToPatch} = 0;%zeros(1,size(qCover{iModeToPatch},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iModeToPatch);

ac_inward = [];
bc_inward = [];
errTrans = [];
newRegArray = [];
bc = [];

idx = 750;  % index of an transition funnel at which to attempt to 'lasso' its cross-section with an inward-facing funnel (i.e. compose the two)

% [ac,existingRegNew,newRegNew] = computePolytopeAtomicController(u0,x0,sys,acNext,reg(aut.q{iModeToPatch}),options.ctrloptions_trans,options.sampSkipFun,xMssExt);
x00 = x0;  u00 = u0;
[ac, c] = computeAtomicController(u00,x00,sys,regMode,ellToCompose,options);
plot(ac.x0,'k',5)
rhoi = double(ac.rho,0);
[res, isectIdx, isectArray] = isinside(ac,reg(aut.q{trans(itrans,1)}),sys);
if ~res
    % NB: the following assumes only one contiguous interval where the funnel left the region.
    % split the funnel into three parts:
    %   - one from the start to the time of the first intersection (acPre)
    %   - another during the interval of intersection (this is created using CBFs)
    %   - another following the intersection (acPost)
    t = ac.x0.getTimeVec();
    acPre = [];
    
    t1 = t(max(isectIdx)+1:end);
    x1 = Traject(t1,double(ac.x0,t1));
    u1 = Traject(t1,double(ac.u0,t1));
    K1 = Traject(t1,double(ac.K,t1));
    P1 = Traject(t1,double(ac.P,t1));
    rho1 = Traject(t1,double(ac.rho,t1));
    Vquad = ac.V(max(isectIdx)+1:end);
    acPost = QuadraticAC(x1,u1,K1,P1,rho1,Vquad,sys);
    
    if false %min(isectIdx) > 1
        t0 = t(1:min(isectIdx)-1);
        x0 = Traject(t0,double(ac.x0,t0));
        u0 = Traject(t0,double(ac.u0,t0));
        
        % Compute a new atomic controller that minimizes the funnel (in contrast to the original method that maximizes it). 
        % The minimization is used here to find a suitable initial condition for the CBF.
        [acPre] = computeAtomicController(u0,x0,sys,regMode,ellToCompose,options,rhoi);
        
    end
    
    % Now, construct the barriers
    t01 = t(min(isectIdx):max(isectIdx));
    x01 = Traject(t01,double(ac.x0,t01));
    u01 = Traject(t01,double(ac.u0,t01));                                        
    
    regMode = reg(aut.q{iModeToPatch});
    
    % define the initial set as the end of the prefix funnel
    if ~isempty(acPre)
        ellToCompose = acPre.ellipsoid(end);
    else
        ellToCompose = acNext.ellipsoid(indexToCompose);
    end

    [ac, bc] = computeConformingFunnel(u01,x01,sys,regMode,ellToCompose,options);

    % Plot stuff
    figure(90)
    axis equal
    hold on
    plot(reg(1),'r')
    plot(reg(2),'g')
    
    plotBarriers(ac,bc)
    
    % pause
    % plot(acNext,sys,90,[],[0,0,1])
    % plot(acPost,sys,90,[],[0,1,0])
    plot(acNext.x0,'k',90)
    plot(ac.x0,'k',90)
    
end
funFail = false;

if ~funFail
    
    plot(ac.x0,'k',3)
    
    % Create a new region based on the last unverified index of the next funnel.
    idxLast = idx + funStepSize;
    
    % add to the set of inward funnels
    ac_inward = [ac_inward; ac];
    
    
    % Create the new region
    newRegVert = [];
    ellAcNext = projection(acNext,sys);
    for j = length(ttmp):-funStepSize:idxLast
        newRegVert = [buildNewRegion(ellAcNext(j), true); newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    newReg = Region(newRegVert(newRegConvHullIdx,:));
    newReg = intersect(newReg,reg(aut.q{iModeToPatch}));
    
    % Subtract the underapproximated reactive funnel
    newRegVert = buildNewRegion(ellAcNext(idxLast), true);
    
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    subReg = Region(newRegVert(newRegConvHullIdx,:));
    subReg = intersect(subReg,reg(aut.q{iModeToPatch}));
    
    newReg = regiondiff(newReg.p,subReg.p);
    
    %reg = [reg; newReg];
    
    % plot it!
    plot(newReg,'m')
    
    newRegVert = extreme(newReg);
    newRegArray = [newRegArray; newReg];
    
    % save the barriers
    bc_inward = [bc_inward; bc];
end
            