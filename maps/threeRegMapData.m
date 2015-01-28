classdef threeRegMapData < MapData
    % Map data for a three-region workspace.
    % One can modify this class to bring a map in from a datafile.
    
    properties(SetAccess = private, GetAccess = public)
        dBloat = 1.2;
        epsilon = 0.1;
        pix2m = 10/(100*2.65);
        
        %LTLMoP-style data sructures
        bndPos = [145, 73]*pix2m;
        bndPoints = [0, 0; 261, 0; 261, 250; 0, 250]*pix2m;
        r1pos = [0, 0]*pix2m;
        r1points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
        r2pos = [175, 100]*pix2m;
        r2points = [0, 0; 90, 0; 90, 134; 0, 134]*pix2m;
        
        reg;
        regBnd;
        regX;
    end
    
    methods
        function obj = threeRegMapData
            obj = obj@MapData;
            
            vReg{1} = [
                obj.r1pos+obj.r1points(1,:);
                obj.r1pos+obj.r1points(2,:);
                obj.r1pos+obj.r1points(3,:);
                obj.r1pos+obj.r1points(4,:)];
            vReg{3} = [
                obj.r2pos+obj.r2points(1,:);
                obj.r2pos+obj.r2points(2,:);
                obj.r2pos+obj.r2points(3,:);
                obj.r2pos+obj.r2points(4,:)];
            vReg{2} = [
                [vReg{3}(2,1) vReg{1}(2,2)];
                vReg{3}(2,:); vReg{3}(1,:);
                vReg{3}(4,:);
                [vReg{1}(1,1) vReg{3}(3,2)];
                vReg{1}(4,:); vReg{1}(3,:);
                vReg{1}(2,:)];
            % vBnd{1} = [vReg{1}(1,:);
            %     vReg{2}(3,:);
            %     vReg{3}(1,:);
            %     vReg{3}(4,:);
            %     vReg{3}(3,:);
            %     vReg{2}(1,:);
            %     vReg{1}(3,:);
            %     vReg{1}(2,:)];  % Bound is based on all the regions
            vBndOuter{1} = [
                vReg{1}(1,:);
                vReg{2}(1,:);
                vReg{3}(3,:);
                vReg{2}(5,:)];
            vBnd = vBndOuter;
            
            % Bloated regions
            % positive bloat: smaller avoid regions
            vRegB{1} = [
                vReg{1}(1,:);
                vReg{1}(2,:)-obj.dBloat*[1 0];
                vReg{1}(3,:)-obj.dBloat*[1 1];
                vReg{1}(4,:)-obj.dBloat*[0 1]];
            vRegB{3} = [
                vReg{3}(1,:)+obj.dBloat*[1 1];
                vReg{3}(2,:)+obj.dBloat*[0 1];
                vReg{3}(3,:);
                vReg{3}(4,:)+obj.dBloat*[1 0]];
            vRegB{2} = [
                [vReg{3}(2,1) vReg{1}(2,2)];
                vReg{3}(2,:)-obj.dBloat*[0 1];
                vReg{3}(1,:)-obj.dBloat*[1 1];
                vReg{3}(4,:)-obj.dBloat*[1 0];
                [vReg{1}(1,1) vReg{3}(3,2)];
                vReg{1}(4,:)+obj.dBloat*[0 1];
                vReg{1}(3,:)+obj.dBloat*[1 1];
                vReg{1}(2,:)+obj.dBloat*[1 0]];
            % negative bloat: larger interior regions
            vRegBN{1} = [
                vReg{1}(1,:);
                vReg{1}(2,:)+obj.dBloat*[1 0];
                vReg{1}(3,:)+obj.dBloat*[1 1];
                vReg{1}(4,:)+obj.dBloat*[0 1]];
            vRegBN{3} = [
                vReg{3}(1,:)-obj.dBloat*[1 1];
                vReg{3}(2,:)-obj.dBloat*[0 1];
                vReg{3}(3,:);
                vReg{3}(4,:)-obj.dBloat*[1 0]];
            vRegBN{2} = [
                [vReg{3}(2,1) vReg{1}(2,2)];
                vReg{3}(2,:)+obj.dBloat*[0 1];
                vReg{3}(1,:)+obj.dBloat*[1 1];
                vReg{3}(4,:)+obj.dBloat*[1 0];
                [vReg{1}(1,1) vReg{3}(3,2)];
                vReg{1}(4,:)-obj.dBloat*[0 1];
                vReg{1}(3,:)-obj.dBloat*[1 1];
                vReg{1}(2,:)-obj.dBloat*[1 0]];
            
            % These regions are defined for checking if funnels are misbehaving... need
            % some conservatism
            vBndOuterB{1} = [
                vReg{1}(1,:)+obj.epsilon*[-1 -1];
                vReg{2}(1,:)+obj.epsilon*[1 -1];
                vReg{3}(3,:)+obj.epsilon*[1 1];
                vReg{2}(5,:)+obj.epsilon*[-1 1]];
            vRegB2{1} = [
                vReg{1}(1,:);
                vReg{1}(2,:)-(obj.dBloat+obj.epsilon)*[1 0];
                vReg{1}(3,:)-(obj.dBloat+obj.epsilon)*[1 1];
                vReg{1}(4,:)-(obj.dBloat+obj.epsilon)*[0 1]];
            vRegB2{3} = [
                vReg{3}(1,:)+(obj.dBloat+obj.epsilon)*[1 1];
                vReg{3}(2,:)+(obj.dBloat+obj.epsilon)*[0 1];
                vReg{3}(3,:);
                vReg{3}(4,:)+(obj.dBloat+obj.epsilon)*[1 0]];
            vRegB2{2} = [
                [vReg{3}(2,1) vReg{1}(2,2)];
                vReg{3}(2,:)-(obj.dBloat+obj.epsilon)*[0 1];
                vReg{3}(1,:)-(obj.dBloat+obj.epsilon)*[1 1];
                vReg{3}(4,:)-(obj.dBloat+obj.epsilon)*[1 0];
                [vReg{1}(1,1) vReg{3}(3,2)];
                vReg{1}(4,:)+(obj.dBloat+obj.epsilon)*[0 1];
                vReg{1}(3,:)+(obj.dBloat+obj.epsilon)*[1 1];
                vReg{1}(2,:)+(obj.dBloat+obj.epsilon)*[1 0]];
            
            vRegBN2{1} = [
                vReg{1}(1,:);
                vReg{1}(2,:)-(obj.epsilon)*[-1 0];
                vReg{1}(3,:)-(obj.epsilon)*[-1 -1];
                vReg{1}(4,:)-(obj.epsilon)*[0 -1]];
            vRegBN2{3} = [
                vReg{3}(1,:)-(obj.epsilon)*[1 1];
                vReg{3}(2,:)-(obj.epsilon)*[0 1];
                vReg{3}(3,:);
                vReg{3}(4,:)-(obj.epsilon)*[1 0]];
            vRegBN2{2} = [
                [vReg{3}(2,1) vReg{1}(2,2)];
                vReg{3}(2,:)-(obj.epsilon)*[0 1];
                vReg{3}(1,:)-(obj.epsilon)*[1 1];
                vReg{3}(4,:)-(obj.epsilon)*[1 0];
                [vReg{1}(1,1) vReg{3}(3,2)];
                vReg{1}(4,:)+(obj.epsilon)*[0 1];
                vReg{1}(3,:)+(obj.epsilon)*[1 1];
                vReg{1}(2,:)+(obj.epsilon)*[1 0]];
            
            % Construct polytopes and hyperplanes
            for i = 1:length(vReg)
                pReg{i} = polytope(vReg{i});
                pRegB{i} = polytope(vRegB{i});
                pRegB2{i} = polytope(vRegB2{i});
                pRegBN{i} = polytope(vRegBN{i});
                pRegBN2{i} = polytope(vRegBN2{i});
                [v,c] = double(pRegBN2{i});
                hRegBN2{i} = hyperplane(v',c');
            end
            
            pBnd{1} = polytope(vBndOuter{1});
            [v,c] = double(pBnd{1});
            hBnd{1} = hyperplane(v',c');
            pBndB{1} = polytope(vBndOuterB{1});
            [v,c] = double(pBndB{1});
            hBndB{1} = hyperplane(v',c');
            
            % Create msspoly arrays
            x = msspoly('x',2);
            for i = 1:length(pReg)
                [H,K] = double(pReg{i});
                mssReg{i} = (H*x(1:2)-K)' + eps*sum(x);
                [H,K] = double(pRegB{i});
                mssRegB{i} = (H*x(1:2)-K)' + eps*sum(x);
                if ~isempty(double(pRegBN{i}))
                    [H,K] = double(pRegBN{i});
                    mssRegBN{i} = (H*x(1:2)-K)' + eps*sum(x);
                else
                    mssRegBN{i} = mssReg{i};
                end
                
                obj.reg{i}.p = pReg{i};
                obj.reg{i}.pB = pRegB{i};
                obj.reg{i}.pB2 = pRegB2{i};
                obj.reg{i}.hBN2 = hRegBN2{i};
                obj.reg{i}.mssExt = -mssReg{i};
                obj.reg{i}.mssExtB = -mssRegBN{i};
                obj.reg{i}.mssInt{1} = [];
                obj.reg{i}.mssIntB{1} = [];
                obj.reg{i}.v = vReg{i};
                obj.reg{i}.vB = vRegB{i};
                obj.reg{i}.vBN = vRegBN{i};
            end
            % TODO: automate the setting of internal regions??
            obj.reg{2}.mssInt{1} = mssReg{1};
            obj.reg{2}.mssInt{2} = mssReg{3};
            obj.reg{2}.mssIntB{1} = mssRegB{1};
            obj.reg{2}.mssIntB{2} = mssRegB{3};
            
            for i = 1:length(pBnd)
                [H,K] = double(pBnd{i});
                mssBnd{i} = (H*x(1:2)-K)' + eps*sum(x);
                
                obj.regBnd{i}.p = pBnd{i};
                obj.regBnd{i}.hB = hBndB{i};
                obj.regBnd{i}.mss = -mssBnd{i};
                obj.regBnd{i}.v = vBnd{i};
            end
            
            % Construct composite region for each transition
            obj.regX{1}.mssExt = obj.regBnd{1}.mss;
            obj.regX{1}.mssExtB = obj.regBnd{1}.mss;
            obj.regX{1}.mssInt{1} = mssReg{3};
            obj.regX{1}.mssIntB{1} = mssRegB{3};
            obj.regX{1}.pExt = obj.regBnd{1}.p;
            obj.regX{1}.hExtB = obj.regBnd{1}.hB;
            obj.regX{1}.pInt(1) = pReg{3};
            obj.regX{1}.pIntB(1) = pRegB{3};
            obj.regX{1}.pIntB2(1) = pRegB2{3};
            
            obj.regX{2} = obj.regX{1};
            obj.regX{2}.mssInt{1} = mssReg{1};
            obj.regX{2}.mssIntB{1} = mssRegB{1};
            obj.regX{2}.pInt(1) = pReg{1};
            obj.regX{2}.pIntB(1) = pRegB{1};
            obj.regX{2}.pIntB2(1) = pRegB2{1};
            
            obj.regX{3} = obj.regX{1};
            obj.regX{4} = obj.regX{2};
            obj.regX{5} = obj.regX{1};
            
            obj.regX{6} = obj.regX{1};
            obj.regX{6}.mssInt{1} = [];
            obj.regX{6}.mssIntB{1} = [];
            obj.regX{6}.pInt(1) = [];
            obj.regX{6}.pIntB(1) = [];
            obj.regX{6}.pIntB2(1) = [];
            
        end
    end
end