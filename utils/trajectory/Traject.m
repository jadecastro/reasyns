
classdef Traject
    
    properties
        pp;
        origbreaks; % TODO: want this to be set when constructing the Traject, and not able to be changed thereafter
    end
    
    methods
        function obj = Traject(arg1,arg2)
            % Constructor
            if isa(arg1,'struct') && isfield(arg1,'form')
                ppobj = arg1;
            elseif isa(arg1,'double') && isa(arg2,'double')
                t = arg1;
                x = arg2;
                if ~all(diff(t) > 0), error('t must be monotonically increasing.'), end
                ppobj = spline(t,x);
            elseif isa(arg1,'double')
                error('representations other than double arrays not implemented.')
            else
                error('input(s) must be either a pp structure or a pair of time and data vectors')
            end
            obj.pp = ppobj;
        end
        
        function [x,t] = double(obj,t)
            %
            if nargin < 2
                t = getTimeVec(obj);
            end
            x = ppval(t,obj.pp);
        end
        
        function [x,t] = downsample(obj,sampSkip)
            %
            [x,t] = double(obj);
            x = downsample(x,sampSkip);
            t = downsample(t,sampSkip);
        end
        
        function [x,t] = downsampleUniformly(obj,sampSkip)
            %
            [x,t] = double(obj);
            [t,indx] = downsampleUniformly(t,sampSkip);
            x = x(:,indx);
        end
        
        function [obj1, obj2] = bisect(obj,Noverlap)
            %
            if nargin < 2
                Noverlap = 0;
            end
            [x,t] = double(obj);
            t1 = t(1:floor(end/2)+Noverlap);
            x1 = x(:,1:floor(end/2)+Noverlap,:);
            t2 = t(floor(end/2):end);
            x2 = x(:,floor(end/2):end,:);
            obj1 = Traject(t1,x1);
            obj2 = Traject(t2,x2);
        end
        
        function obj = horzcat(obj1,obj2)
            % 
            t1 = getTimeVec(obj1);
            t2 = getTimeVec(obj2);
            if (t1(1) < t2(1)) || (t1(end) > t2(end)), 
                warning('Time vector of first argument not within the range of that of the second argument. You may get large extrapolation errors!')
            end
            X = ppval(t1,obj1.pp);
            Y = ppval(t1,obj2.pp);
            ndx = ndims(X);  ndy = ndims(Y);
            
            if ndx ~= ndy, error('Dimensions must match!'), end
            S1.type = '()';
            for i = 1:ndx-1
                S1.subs{i} = 1:size(X,i);
            end
            S2 = S1;
            S2.subs{2} = 1:(size(X,2)+size(Y,2));
            Z = [];
            for i = 1:length(t1)
                S1.subs{ndx} = i;  
                S2.subs{ndx} = i;
                tmpX = subsref(X,S1);  tmpY = subsref(Y,S1);
                tmpZ = [tmpX tmpY];
                Z = subsasgn(Z,S2,tmpZ);
            end
            obj = Traject(t1,Z);
        end
        
        function obj = vertcat(obj1,obj2)
            %
            t1 = getTimeVec(obj1);
            t2 = getTimeVec(obj2);
            if (t1(1) < t2(1)) || (t1(end) > t2(end)),
                warning('Time vector of first argument not within the range of that of the second argument. You may get large extrapolation errors!')
            end
            X = ppval(t1,obj1.pp);
            Y = ppval(t1,obj2.pp);
            ndx = ndims(X);  ndy = ndims(Y);
            
            if ndx ~= ndy, error('Dimensions must match!'), end
            S1.type = '()';
            for i = 1:ndx-1
                S1.subs{i} = 1:size(X,i);  
            end
            S2 = S1;
            S2.subs{1} = 1:(size(X,1)+size(Y,1));
            Z = [];
            for i = 1:length(t1)
                S1.subs{ndx} = i;  
                S2.subs{ndx} = i;
                tmpX = subsref(X,S1);  tmpY = subsref(Y,S1);
                tmpZ = [tmpX; tmpY];
                Z = subsasgn(Z,S2,tmpZ);
            end
            obj = Traject(t1,Z);
        end
        
        function obj = timecat(obj1,obj2)
            %
            [x1,t1] = double(obj1);
            [x2,t2] = double(obj2);
            if t2(1) < t1(end), error('Resulting time vector must be monotonic when the two input vectors are concatenated.'), end
            ndx1 = ndims(x1); ndx2 = ndims(x2);
            if ndx1 ~= ndx2, error('Dimensions must match!'), end
            for i = 1:ndx1-1
                if size(x1,i) ~= size(x2,i), error('size of data arrays must match.'), end
            end
            t = [t1 t2];
            x = cat(ndx1,x1,x2);
            [t,indxkeep] = unique(t);
            
            %TODO: verify that data for each overlapping element are indeed
            %the same
            if length(indxkeep) < size(x,ndx1)
                S1.type = '()';
                for i = 1:ndx1-1
                    S1.subs{i} = 1:size(x1,i);
                end
                S2 = S1;
                z = [];
                for i = 1:length(indxkeep)
                    S1.subs{ndx1} = indxkeep(i);
                    S2.subs{ndx1} = i;
                    tmpz = subsref(x,S1);
                    z = subsasgn(z,S2,tmpz);
                end
            else
                z = x;
            end
            obj = Traject(t,z);    
        end
        
        function obj = rdivide(obj1,obj2)
            %
            if obj2.pp.dim > 1, error('divisor (2nd arg) dimension must be one'); end
            t1 = getTimeVec(obj1);
            t2 = getTimeVec(obj2);
            if (t1(1) < t2(1)) || (t1(end) > t2(end)),
                warning('Time vector of first argument not within the range of that of the second argument. You may get large extrapolation errors!')
            end
            X = ppval(t1,obj1);
            Y = ppval(t1,obj2);
            Z = X./Y;
            obj = spline(t1,Z);
        end
        
        function obj = times(obj1,obj2)
            %
            if obj2.pp.dim > 1, error('multiplier (2nd arg) dimension must be one'); end
            t = getTimeVec(obj1);
            X = ppval(t,obj1);
            Y = ppval(t,obj2);
            Z = X.*Y;
            obj = spline(t,Z);
        end
        
        function t = getTimeVec(obj)
            %
            t = obj.pp.breaks;
        end
        
        function len = length(obj)
            %
            len = length(obj.pp.breaks);
        end
        
        function plot(obj,color,fignum)
            %
            if nargin == 3 && isempty(color)
                color = 'k';
            end
            if nargin < 3
                fignum = 1;
            end
            if nargin < 2
                color = 'k';
            end
            [x0,t] = double(obj);
            % TODO: plot using coord transformation matrix
            figure(fignum)
             plot(x0(1,:),x0(2,:),'LineWidth',2)
            drawnow
        end
    end
    
end