function autNew = renumberStates(aut)

currentStates = unique([aut.trans{:}]);

for idx = 1:length(currentStates)
    autNew.q{idx} = aut.q{currentStates(idx)};
end

count = 0;
for idx = 1:length(aut.trans)
    if ~isempty(aut.trans{idx})
        count = count+1;
        autNew.trans{count} = [find(currentStates == aut.trans{idx}(1)) find(currentStates == aut.trans{idx}(2))];
    end
end
    