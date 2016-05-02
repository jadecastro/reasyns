function [indexTransVect] = findTransitionsWithNonRepeatingRegions(aut,trialState)

trans = vertcat(aut.trans{:});

indexTransVectFull = find(trans(:,1)==trialState)';
if ~isempty(indexTransVectFull)
    [~,indexPostVect] = intersect(vertcat(aut.state{:}),trans(indexTransVectFull,2));
    [~,indexOfUniqueIndexVect] = unique(vertcat(aut.label{indexPostVect}));  % eliminate any post states that are redundant wrt region labels
    indexUniquePostVect = indexPostVect(indexOfUniqueIndexVect);
    [~,indexOfPostIndexVect] = setdiff(vertcat(aut.label{indexUniquePostVect}),aut.label{vertcat(aut.state{:}) == trialState}); % sliminate any transitions leading to the same region
    [~,indexOfTransVect] = intersect(trans(indexTransVectFull,2),vertcat(aut.state{indexUniquePostVect(indexOfPostIndexVect)}));
    indexTransVect = indexTransVectFull(indexOfTransVect);
else
    indexTransVect = [];
end
%indexPostVect = trans(indexTransVect,2);

% skip if any atomic controllers have already been created