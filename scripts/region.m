
classdef region < handle
    
    properties (SetObservable)
        numRegions;
        
        % region data
        v;
        p;
        mssExt;
        mssInt;
        mssBnd;
        
        p_bnd;
        mss_bnd;
    end
    
    properties (SetAccess = private)
        epsilon;
        inflated = false;
        parents;
        children;
    end

    methods
        function obj = region(v,epsilon,parentreg)
            % Constructor
            % Note: vertices must be specified in counterclockwise order
            
            if ~nargin 
                return
            end
            
            if isnumeric(v)
                if nargin > 2
                    v = parentreg.v;
                end
                obj.v = v;
            else
                error('input arg must be an n-by-2 array of vertices.')
            end
            
            configureRegion(obj,v);
            
            if nargin > 1
                obj.epsilon = epsilon;
                inflate(obj);
            end
            addlistener(obj,'v','PostSet',@region.handleSelfVertexChanges);
            if nargin > 2
                % TODO: not very robust-- set-up an event trigger instead?
                %                 addlistener(parentreg,'p','PostSet',@(src,evt)handleOtherVertexChanges(this,src,evt));
                addlistener(parentreg,'p','PostSet',@obj.handleOtherVertexChanges);
            end
        end
        
        function setNumRegions(regobj,val)
            regobj.numRegions = val;
        end

        function trueRep(regBnd,regArray)
            % regRes = trueRep(regBnd,regArray)
            % Compute an array of convex polygons to represent the true region
            % regRes such that regRes = regBnd \ regArray.
            pArray = [];
            for i = 1:length(regArray)
                pArray = [pArray; regArray(i)];
            end
            regBnd.p = regiondiff(regBnd.p,pArray);
        end
        
        function inflate(regobj)
            % regRes = inflate(regobj,epsilon)
            % Inflate the regobj by a value of epsilon (negative values of
            % epsilon result in a deflation of regobj)
            if regobj.epsilon ~= 0 && ~regobj.inflated
                x = msspoly('x',2);
                [H,K] = double(regobj.p);
                if length(regobj.p) > 1
                    for i = 1:length(regobj.p)
                        tmp = polytope(H{i},K{i}+regobj.epsilon);
                        if ~isempty(double(tmp))
                            regobj.p(i) = tmp;
                        end
                        mssReg = (H{i}*x(1:2)-(K{i}+regobj.epsilon))' + eps*sum(x);
                        regobj.mssExt{i} = -mssReg;
                    end
                else
                    regobj.p = polytope(H,K+regobj.epsilon);
                    mssReg = (H*x(1:2)-(K+regobj.epsilon))' + eps*sum(x);
                    regobj.mssExt = -mssReg;
                end
                
                % deflate the boundary 
                [H,K] = double(regobj.p_bnd);
                regobj.p_bnd = polytope(H,K-regobj.epsilon);
                mssReg = (H*x(1:2)-(K-regobj.epsilon))' + eps*sum(x);
                regobj.mss_bnd = -mssReg;
                
                %                 regobj.v = extreme(regobj.p);
                %                 regobj.v = [];
                regobj.inflated = true; % flag as inflated
            end
        end
        
        function h = gethyperplane(regobj)
            %
            [v,c] = double(regobj.p);
            h = hyperplane(v',c');
        end
        
        function mss = getMSSpoly(regobj)
            %
            
        end
        
        function reg = intersect(reg1,reg2)
            %  find the intersection of two region objects, returning a new region object
            
            p = intersect(reg1.p,reg2.p);
            verts = extractOrderedVertsFromPolytope(p);
            reg = region(verts);
        end
        
        function reg = union(reg1,reg2)
            %  find the union of two region objects, returning a new region object
            
            p = union([reg1.p,reg2.p]);
            verts = extractOrderedVertsFromPolytope(p);
            reg = region(verts);
        end
        
        function isect = isinside(regobj,sys,q,sampSkip)
            % check whether or not any of the states values supplied in q do not lie in a given region.
            
            H = sys.params.H;
            if isa(q,'traject')
                q = double(q);
            end
            
            if nargin > 3
                q = downsample(q,sampSkip);
            end
            
            if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
            
            for i = 1:size(q,1)
%                 try
                isect = isinside(vertcat(regobj.p),H*q(i,1:length(H))');
%                 catch
%                     keyboard
%                 end
                if ~isect
                    return
                end
            end
        end
        
        function regDiff = regdiff(regobj,regArray)
            % 
            % NB: assumes regobj is convex.
            polyArr = [];
            for i = 1:length(regArray)
                for j = 1:length(regArray(i).p)
                    polyArr = [polyArr; regArray(i).p(j)];
                end
            end
            pDiff = regiondiff(regobj.p,polyArr);
            vDiff = [];
            for i = 1:length(pDiff)
                vDiffNew = extreme(pDiff(i));
                if i > 1
                    if verLessThan('matlab','7.15')
                        P1.x = vDiff(:,1); P1.y = vDiff(:,2);
                        P2.x = vDiffNew(:,1); P2.y = vDiffNew(:,2);
                        P3 = PolygonClip(P1,P2,3);
                        a = P3.x; b = P3.y;
                    else
                        [a,b] = polybool('union',vDiff(:,1),vDiff(:,2),vDiffNew(:,1),vDiffNew(:,2));
                    end
                    vDiff = [a b];
                else
                    vDiff = vDiffNew;
                end
            end
            if verLessThan('matlab','7.15') % use the DIY version
                a = vDiff(:,1); b = vDiff(:,2);
                amean=mean(vDiff(:,1));  bmean=mean(vDiff(:,2));
                [~, indx] = sort(atan2(a-amean,b-bmean));
                a = a(indx);  b = b(indx);
            else
                [a,b] = poly2ccw(vDiff(:,1),vDiff(:,2));
            end
                
            regDiff = region([a b]);
        end
        
        function [regSafe] = getReg(regobj,regbnd,aut,imode)
            % Safe region for the current imode
            count = 0;
            for j = 1:length(regobj)
                if aut.q{imode} ~= j
                    count = count+1;
                    regAvoid(count) = regobj(j);
                end
            end
            regSafe = regdiff(regbnd,regAvoid);
        end
        
        function [regSafe] = getRegTrans(regobj,regbnd,aut,itrans)
            % Safe regions for the given itrans
            count = 0;
            for j = 1:length(regobj)
                if (aut.q{aut.trans{itrans}(1)} ~= j) && (aut.q{aut.trans{itrans}(2)} ~= j)
                    count = count+1;
                    regAvoid(count) = regobj(j);
                end
            end
            regSafe = regdiff(regbnd,regAvoid);
        end
        
        function plot(regobj,varargin)
            %
            plot(regobj.p,varargin{:});
        end
        
%         function subsasgn
%             %
%             error('wip')
%         end
  
        function handleOtherVertexChanges(regobj,src,evt)
            %
            disp(['reconfiguring the child region'])
            epsilon = regobj.epsilon;
            regobj.v = evt.AffectedObject.v;
            regobj.p = evt.AffectedObject.p;
            regobj.epsilon = epsilon;
            regobj.inflated = false;
%             configureRegion(regobj,evt.AffectedObject.v);
            inflate(regobj);
        end

    end
        
    methods (Static)
        
        function handleSelfVertexChanges(src,evt)
            %
            disp(['configuring the current region'])
            configureRegion(evt.AffectedObject,evt.AffectedObject.v);
        end            
    end
    
    methods (Access = private)
        
        function configureRegion(regobj,v)
            regobj.p = polytope(regobj.v);
            %             regobj.v = extreme(regobj.p);
            
            % polytope returns only a single convex region; we need to find an array
            % of non-convex regions
            numOrigVerts = size(v,1);
            indxNC = [];
            h = gethyperplane(regobj);
            for i = 1:numOrigVerts
                if ~any(intersect(ellipsoid(v(i,:)',eye(2)*1e-20),h)), indxNC = [indxNC; i]; end
            end
            p1 = [];
            for i = 1:numOrigVerts
                clear v1
                if ~ismember(i,indxNC) && ismember(mod(i,numOrigVerts)+1,indxNC)
                    v1(1,:) = v(i,:);
                    v1(2,:) = v(mod(i,numOrigVerts)+1,:);
                    for j = 1:numOrigVerts
                        v1(2+j,:) = v(mod(i+j,numOrigVerts)+1,:);
                        if ~ismember(mod(i+j,numOrigVerts)+1,indxNC), break; end
                    end
                    p1 = [p1; polytope(v1)];
                end
            end
            if ~isempty(p1)
                trueRep(regobj,p1);
            end
            
            % compute an mss object for the polytope
            [H,K] = double(regobj.p);
            x = msspoly('x',2);
            if length(regobj.p) > 1
                for i = 1:length(regobj.p)
                    mssReg = (H{i}*x(1:2)-K{i})' + eps*sum(x);
                    regobj.mssExt{i} = -mssReg;
                end
            else
                mssReg = (H*x(1:2)-K)' + eps*sum(x);
                regobj.mssExt = -mssReg;
            end
            
            % form a convex hull of the polytope
            bnd = hull(regobj.p);
            regobj.p_bnd = bnd;
            
            % compute an mss object for the convex hull
            [H,K] = double(regobj.p_bnd);
            mssReg = (H*x(1:2)-K)' + eps*sum(x);
            regobj.mssBnd = -mssReg;

        end
        
    end
end