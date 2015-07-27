
function [aut, transTmp, transTmpNew, transTmpNewWithoutSelfLoops, transWithoutSelfLoops] = processAutFileFastSlow(autfname)
% Return an automaton structure from data that is read in from an aut file.

alphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

% read from the .aut file
% fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/complex_map/complex_map.aut','r');
% fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/light_painting_letters/light_painting.aut','r');
fid = fopen(autfname,'r');

idx = 1;
clear A
tline = fgetl(fid);
A{idx} = tline;
while ischar(tline)
    idx = idx+1;
    tline = fgetl(fid);
    A{idx} = tline;
end
fclose(fid);

regPairsNonUnique = [];
currTmp = [];
nextTmp = [];
rankTmp = [];
nextTmpIdx = 0;
rankAndNext = [];
savedCurr = [];
testCurrState = [];
for idx = 1:length(A)
    Key   = 'State';
    Index = strfind(A{idx}, Key);
    if ~isempty(Index)
        stateVal = sscanf(A{idx}(Index(1) + length(Key):end), '%g', 1);
        testCurrState = [testCurrState; stateVal];
        
        Key   = ['with rank'];
        Index = strfind(A{idx}, Key);
        rankVal = sscanf(A{idx}(Index(1) + length(Key):end), '%g', 1);
        
        currTmp = [currTmp; stateVal];  % columns correspond to lines in the aut file. Originally, these are the same as the state designations, but will be re-mapped when we flatten all of the 'duplicate' states.
        rankTmp = [rankTmp; rankVal];

        
        for currReg = 1:length(alphabet)
            Key   = ['rob1_',alphabet{currReg},'_rc:'];
            Index = strfind(A{idx}, Key);
            if sscanf(A{idx}(Index(1) + length(Key):end), '%g', 1)
                break
            end
        end
        for nextReg = 1:length(alphabet)
            Key   = ['rob1_',alphabet{nextReg},':'];
            Index = strfind(A{idx}, Key);
            if sscanf(A{idx}(Index(1) + length(Key):end), '%g', 1)
                break
            end
        end
    end
    
    % store the set of successors to the current state (alternates with the current state-populating 'if' statement above)  
    Key   = 'With successors :';
    Index1 = strfind(A{idx}, Key);
    if ~isempty(Index1)
        succVal = sscanf(A{idx}(Index1(1) + length(Key):end), '%g', 1);
        indx = 1;
        for k = strfind(A{idx}, ',')
            indx = indx+1;
            succVal = [succVal sscanf(A{idx}(k+1:end), '%g', 1)];
        end
        nextTmpIdx = nextTmpIdx+1;
        nextTmp{nextTmpIdx} = succVal;
    
        % if a state in the list has the same region label with the same rank as an existing state, map this new state to the existing one.
        %         succVal
        if ~isempty(succVal)
            [mbrTrue,idxRank] = ismember([currReg rankVal],rankAndNext,'rows');
            if mbrTrue
                currTmp(end) = savedCurr(idxRank);
            else
                rankAndNext = [rankAndNext; currReg rankVal];
                savedCurr = [savedCurr; stateVal];
            end
            if currReg == 1
%                 [currReg nextReg]
            end
        else
            currTmp(end) = NaN;  % mark states where the system can force a violation of the environment safety conditions as 'NaN'   
        end
        regPairsNonUnique = [regPairsNonUnique; currReg, nextReg];
    end
end
if ~all(testCurrState(2:end) - testCurrState(1:end-1))
    error('the state vector must be sorted and always increasing!')
end

% map the list of next states to the reduced set and construct the aut structure
clear aut
count = 0;
transTmp = [];
for idx = 1:length(currTmp)
% for idx = unique(currTmp)'+1
    aut.q{idx} = regPairsNonUnique(idx,1);
    tmpnxt = [];
    for j = nextTmp{idx}
        tmpnxt = [tmpnxt; currTmp(j+1)];
        %             j
        %             nextTmp{idx}(j)
    end
    nextTmp{idx} = tmpnxt';
    nextTmp{idx} = unique(nextTmp{idx});
    for j = nextTmp{idx}
        if ~isnan(j)
            count = count+1;
            aut.trans{count} = [(currTmp(idx)+1) (j+1)];
            transTmp = [transTmp; regPairsNonUnique(currTmp(idx)+1,1) regPairsNonUnique(j+1,1)];
        end
        if regPairsNonUnique(currTmp(idx)+1,1) == 1 && ~isnan(j)
            [regPairsNonUnique(currTmp(idx)+1,1) regPairsNonUnique(j+1,1)]  % display some of the results, in terms of the regions, as a sanity check
        end
    end
end

% make everything unique
tmp = [];
for idx = 1:length(aut.trans)
    tmp = [tmp; aut.trans{idx}];
end
tmp = unique(tmp,'rows');
aut.trans = [];
transWithoutSelfLoops = [];
transTmpNew = [];
transTmpNewWithoutSelfLoops = [];
for idx = 1:size(tmp,1)
    aut.trans{idx} = tmp(idx,:);
    transTmpNew = [transTmpNew; regPairsNonUnique(aut.trans{idx},1)'];
    
    % modify the transitions and state designations removing duplicate region transitions, regardless of rank  
%     regPairsNonUnique(transWithoutSelfLoops(regPairsNonUnique(transWithoutSelfLoops(:,2))==1,1))

    
    % create a preview of the transitions without self-loops
    if aut.q{aut.trans{idx}(1)} ~= aut.q{aut.trans{idx}(2)}
        transWithoutSelfLoops = [transWithoutSelfLoops; aut.trans{idx}];
        transTmpNewWithoutSelfLoops = [transTmpNewWithoutSelfLoops; regPairsNonUnique(aut.trans{idx},1)'];
    end
end

% get rid of self-loops and flatten any of those that are somewhere followed by an outgoing transition from that state. Note: tmp and aut.trans have the same nubmer of entries. 
tmp = [];
for idx = 1:length(aut.trans)
    done = false;
    newTrans = aut.trans{idx};
    if ~isempty(aut.trans{idx})
        if (aut.q{aut.trans{idx}(1)} == aut.q{aut.trans{idx}(2)})
            if (aut.trans{idx}(1) ~= aut.trans{idx}(2))
                disp('found...')
                [aut.q{aut.trans{idx}(1)} aut.q{aut.trans{idx}(2)}]
                [aut.trans{idx}(1) aut.trans{idx}(2)]
                transPre = aut.trans{idx}(1);
                transPost = aut.trans{idx}(2);
                transPostVec = [];
                tmpSav = [];
                while ~done
                    [done, transPost, transPostLoop] = findOutgoingPost(transPost,aut);
                    transPostVec = [transPostVec; transPost];
                    transPost = transPostLoop;
                    if ismember(transPostLoop, tmpSav)
                        break
                    end
                    tmpSav = [tmpSav; transPostLoop];
                end
                newTrans = [transPre*ones(length(transPostVec),1) transPostVec]
            else
                newTrans = [];
            end
        end
    end
    tmp = [tmp; newTrans];
end
tmp = unique(tmp,'rows');
tmp = tmp(ismember(tmp(:,1),tmp(:,2)),:);
aut.trans = [];
for idx = 1:length(tmp)
    aut.trans{idx} = tmp(idx,:);
end

function [done, transPost, transPost2] = findOutgoingPost(transPost,aut)

done = true;
newTransPost = [];
transPost2 = [];
for idx = 1:length(aut.trans)
    if ~isempty(aut.trans{idx})
        if ((aut.trans{idx}(1) == transPost) && (aut.trans{idx}(2) ~= transPost))
            transPost1 = aut.trans{idx}(2);
            if (aut.q{aut.trans{idx}(1)} == aut.q{aut.trans{idx}(2)})
                disp('1...')
                [aut.q{aut.trans{idx}(1)} aut.q{aut.trans{idx}(2)}]
                [aut.trans{idx}(1) aut.trans{idx}(2)]
                transPost2 = transPost1;
                done = false;
            else
                disp('2...')
                [aut.q{aut.trans{idx}(1)} aut.q{aut.trans{idx}(2)}]
                [aut.trans{idx}(1) aut.trans{idx}(2)]
                newTransPost = [newTransPost; transPost1];
            end
        end
    end
end
if ~isempty(newTransPost)
    transPost = newTransPost;
end


