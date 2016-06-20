
classdef Region < handle
    
    properties (SetObservable)
        numRegions;
        
        % region data
        name;
        v;
        p;
        
        hpp;
        
        mssExt;
        mssInt;
        mssBnd;
        calibMatrix;  % either use a calibration matrix or a simple scalar (of units pixels per meter) to map the raw region data into lab units
        
        p_bnd;
        mss_bnd;
    end
    
    properties (SetAccess = private)
        epsilon;
        isInflated = false;
        isScaled = false;
        parents;
        children;
    end

    methods
        function obj = Region(name,v,calibMatrixOrScalar,epsilon,parentreg)
            % Constructor
            % Note: vertices must be specified in counterclockwise order
            
            if ~nargin 
                return
            end
            
            obj.name = name;
            
            % The calibration matrix assumes an SE(2) system
            if nargin > 2 && ~isempty(calibMatrixOrScalar)
                if isscalar(calibMatrixOrScalar)
                    obj.calibMatrix = diag([1/calibMatrixOrScalar*[1 1] 1]);
                end
                obj.calibMatrix = calibMatrixOrScalar;
                    
            else
                obj.calibMatrix = 1;
            end
            
            if isnumeric(v)
                if nargin > 4
                    v = parentreg.v;
                end
                obj.v = v;
            else
                error('input arg must be an n-by-2 array of vertices.')
            end
                      
            scaleRegion(obj);
            
            configureRegion(obj,v);
            
            if nargin > 3
                obj.epsilon = epsilon;
                inflate(obj);
            end
            addlistener(obj,'v','PostSet',@Region.handleSelfVertexChanges);
            if nargin > 4
                % TODO: not very robust-- set-up an event trigger instead?
                %                 addlistener(parentreg,'p','PostSet',@(src,evt)handleOtherVertexChanges(this,src,evt));
                addlistener(parentreg,'p','PostSet',@obj.handleOtherVertexChanges);
            end
        end
        
        function setNumRegions(obj,val)
            obj.numRegions = val;
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
        
        function inflate(obj)
            % regRes = inflate(obj,epsilon)
            % Inflate the region by a value of epsilon (negative values of epsilon result in a deflation of region)
            if obj.epsilon ~= 0 && ~obj.isInflated
                x = msspoly('x',2);
                [H,K] = double(obj.p);
                if length(obj.p) > 1
                    for i = 1:length(obj.p)
                        tmp = polytope(H{i},K{i}+obj.epsilon);
                        if ~isempty(double(tmp))
                            obj.p(i) = tmp;
                        end
                        mssReg = (H{i}*x(1:2)-(K{i}+obj.epsilon))' + eps*sum(x);
                        obj.mssExt{i} = -mssReg;
                    end
                else
                    obj.p = polytope(H,K+obj.epsilon);
                    mssReg = (H*x(1:2)-(K+obj.epsilon))' + eps*sum(x);
                    obj.mssExt = -mssReg;
                end
                
                % deflate the boundary 
                [H,K] = double(obj.p_bnd);
                obj.p_bnd = polytope(H,K-obj.epsilon);
                mssReg = (H*x(1:2)-(K-obj.epsilon))' + eps*sum(x);
                obj.mss_bnd = -mssReg;
                
                obj.v = extreme(obj.p);
                %                 obj.v = [];
                obj.isInflated = true; % flag as inflated
            end
        end
        
        function h = gethyperplane(obj)
            %
            [v,c] = double(obj.p);
            h = hyperplane(v',c');
        end
        
        function mss = getMSSpoly(obj)
            %
            
        end
        
        function reg = intersect(reg1,reg2)
            %  find the intersection of two region objects, returning a new region object
            
            p = intersect(reg1.p,reg2.p);
            verts = extractOrderedVertsFromPolytope(p);
            if ~isempty(verts)
                reg = Region([reg1.name,'_intersect_',reg2.name],verts);
            else
                reg = [];
            end
        end
        
        function reg = union(reg1,reg2)
            %  find the union of two region objects, returning a new region object
            
            p = union([reg1.p,reg2.p]);
            verts = extractOrderedVertsFromPolytope(p);
            
            reg = Region([reg1.name,'_union_',reg2.name],verts);
            
            % prevent re-scaling the scaled vertices
            reg.calibMatrix = reg1.calibMatrix;
            reg.isScaled = reg1.isScaled;
        end
        
        function reg = regiondiff(reg1,reg2)
            %  find the difference between two region objects, returning a new region object
            
            p = regiondiff(reg1.p,reg2.p);
            verts = extractOrderedVertsFromPolytope(p);
            reg = Region([reg1.name,'_without_',reg2.name],verts,reg1.calibMatrix);
        end
        
        function isect = isinside(obj,sys,q,sampSkip)
            % check whether or not any of the states values supplied in q do not lie in a given region.
            
            if size(q,2) == 1
                q = q';
            end
            
            if isa(q,'Traject')
                q = double(q);
            end
            
            if nargin > 3
                q = downsample(q,sampSkip);
            end
            
            %if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
            
            for i = 1:size(q,1)
%                 try
                output = sys.state2SEconfig([],q(i,:),[]);
                isect = isinside(vertcat(obj.p),output(1:2));
%                 catch
%                     keyboard
%                 end
                if ~isect
                    return
                end
            end
        end
        
        function [res, xFail] = regionContainsEllipsoid(obj,sys,ell)
            % Checks whether or not the funnel contains a given ellipsoid.
            % Returns 'true' if contained and 'false' otherwise.
            
            xFail = [];
            
            ellTestProj = obj.projection(sys,ell);  % Returns empty if not an orthogonal basis (often the case if nonlinear)
            
            if ~isempty(ellTestProj) && sys.isOutputLinear
                res = ~any(intersect(ellTestProj,obj.hpp,'u'));
                
            else  % use a sample-based approximation with no guarantee that a declared containment is actually true
                d = 100;
                [c,Q] = double(ell);
                x = sampleEllipsoidBoundary(c,Q,d);
                
                res = true;
                
                % Loop over all randomly-selected points
                for i = 1:size(x,1)
                    if ~obj.isinside(sys,x(i,:))
                        res = false;
                        xFail = [xFail; x(i,:)];
                        if nargout < 2
                            return;
                        end
                    end
                end
            end
        end
        
        function ellProj = projection(obj,sys,ell)
            
            ellProj = [];
            
            [X,~] = double(ell);            
            [~,~,H] = getRegNonRegStates(sys,[],X,[]);

            try
                ellProj = projection(ell,H(1:2,:)');
            catch ME
                if ~isempty(strfind(ME.message,'PROJECTION'))  % The state-to-output mapping is non-ortogonal; return empty
                    return
                end
            end
        end
        
        function regDiff = regdiff(obj,regArray)
            % 
            % NB: assumes the region is convex.
            polyArr = [];
            for i = 1:length(regArray)
                for j = 1:length(regArray(i).p)
                    polyArr = [polyArr; regArray(i).p(j)];
                end
            end
            pDiff = regiondiff(obj.p,polyArr);
            vDiff = [];
            for i = 1:length(pDiff)
                vDiffNew = extreme(pDiff(i));
                if i > 1
                    if true %verLessThan('matlab','7.15')
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
            if true %verLessThan('matlab','7.15') % use the DIY version
                a = vDiff(:,1); b = vDiff(:,2);
                amean=mean(vDiff(:,1));  bmean=mean(vDiff(:,2));
                [~, indx] = sort(atan2(a-amean,b-bmean));
                a = a(indx);  b = b(indx);
            else
                [a,b] = poly2ccw(vDiff(:,1),vDiff(:,2));
            end
            
            regDiff = Region('Diffed_',[a b],obj.calibMatrix);
        end
        
        function [regSafe] = getReg(obj,regbnd,aut,imode)
            % Safe region for the current imode
            count = 0;
            for j = 1:length(obj)
                if aut.label{vertcat(aut.state{:}) == imode} ~= j
                    count = count+1;
                    regAvoid(count) = obj(j);
                end
            end
            regSafe = regdiff(regbnd,regAvoid);
        end
        
        function [regSafe] = getRegTrans(obj,regbnd,aut,itrans)
            % Safe regions for the given itrans
            count = 0;
            for j = 1:length(obj)
                if (aut.label{vertcat(aut.state{:}) == aut.trans{itrans}(1)} ~= j) && (aut.label{vertcat(aut.state{:}) == aut.trans{itrans}(2)} ~= j)
                    count = count+1;
                    regAvoid(count) = obj(j);
                end
            end
            regSafe = regdiff(regbnd,regAvoid);
        end
        
        function plot(obj,varargin)
            %
            plot(obj.p,varargin{:});
        end
        
%         function subsasgn
%             %
%             error('wip')
%         end
  
        function handleOtherVertexChanges(obj,src,evt)
            %
            disp(['reconfiguring the child region'])
            epsilon = obj.epsilon;
            obj.v = evt.AffectedObject.v;
            obj.p = evt.AffectedObject.p;
            obj.epsilon = epsilon;
            obj.isInflated = false;
%             configureRegion(obj,evt.AffectedObject.v);
            inflate(obj);
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
        
        function configureRegion(obj,v)
            obj.p = polytope(obj.v);
            %             obj.v = extreme(obj.p);
            
            % polytope returns only a single convex region; we need to find an array
            % of non-convex regions
            %TODO: fix
            if false
                numOrigVerts = size(v,1);
                indxNC = [];
                obj.hpp = gethyperplane(obj);
                for i = 1:numOrigVerts
                    if ~any(intersect(ellipsoid(v(i,:)',eye(2)*1e-20),obj.hpp)), indxNC = [indxNC; i]; end
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
                    trueRep(obj,p1);
                end
            end
            
            % compute an mss object for the polytope
            [H,K] = double(obj.p);
            x = msspoly('x',2);
            if length(obj.p) > 1
                for i = 1:length(obj.p)
                    mssReg = (H{i}*x(1:2)-K{i})' + eps*sum(x);
                    obj.mssExt{i} = -mssReg;
                end
            else
                mssReg = (H*x(1:2)-K)' + eps*sum(x);
                obj.mssExt = -mssReg;
            end
            
            % form a convex hull of the polytope
            bnd = hull(obj.p);
            obj.p_bnd = bnd;
            
            % compute an mss object for the convex hull
            [H,K] = double(obj.p_bnd);
            mssReg = (H*x(1:2)-K)' + eps*sum(x);
            obj.mssBnd = -mssReg;

        end
        
        function scaleRegion(obj)
            % scale the region vertices appropriately.

            if ~obj.isScaled
                vNew = [];
                for j = 1:size(obj.v,1)
                    vtmp = inv(obj.calibMatrix)*[obj.v(j,:)';1];
                    vNew = [vNew; vtmp(1:2)'];
                end
                obj.v = vNew;
            end

            obj.isScaled = true;
        end
        
    end
end