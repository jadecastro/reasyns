function qCover = getCoverPts(vReg,trans,itrans,Ncover,H,n,limsNonRegState)

if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
if n == length(H), error('Empty non-region state vector not yet supported. Sorry.'); end
if size(limsNonRegState,2) ~= n - length(H), error('Dimension mismatch in limsNonRegState.'); end

% rand('seed',0);

minNonRegStates = limsNonRegState(1,:);
maxNonRegStates = limsNonRegState(2,:);

for i = 1:Ncover
    
    randRegStates = rand(1,length(H));
    randNonRegStates = rand(1,n-length(H));
    tmpRegStates = (H\(min(vReg{trans{itrans}(1)}) + (max(vReg{trans{itrans}(1)}) - min(vReg{trans{itrans}(1)})).*randRegStates*H)')';
    tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*randNonRegStates;
    
    while ~inpolygon(tmpRegStates(1),tmpRegStates(2),vReg{trans{itrans}(1)}(:,1),vReg{trans{itrans}(1)}(:,2))
        randRegStates = rand(1,length(H));
        randNonRegStates = rand(1,n-length(H));
        tmpRegStates = (H\(min(vReg{trans{itrans}(1)}) + (max(vReg{trans{itrans}(1)}) - min(vReg{trans{itrans}(1)})).*randRegStates*H)')';
        tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*randNonRegStates;
    end
    qCover(i,:) = [tmpRegStates tmpNonRegStates];
end

