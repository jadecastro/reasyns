
pix2m = 1;
% calibMatrix =   [ 33.50049738,  -0.89279945,  66.48435372,
%                   -0.27159749, -30.69295068,  67.89619503,
%                    0.        ,   0.        ,   1.        ];

% update the existing region, add region to region list, and assign it a new state
% reg(aut.q{iModeToPatch}) = existingReg;
reg(length(reg)+1) = newRegVert;
newReg = region(newRegVert);
aut.q{length(aut.q)+1} = max([aut.q{:}]) + 1;

vertsToWrite = 1/pix2m*calibMatrix*[newReg.v ones(size(newReg.v,1),1)]';
vertsToWrite(3,:) = [];
vertsToWrite = vertsToWrite';

% transform to map coordinates
% for idx = 1:size(vertsToWrite,1)
%     tmpVert(idx,:) = 1/pix2m*(calibMatrix*[vertsToWrite(idx,1:2)'; 1])';
% end
% vertsToWrite = tmpVert(:,1:2);

% compose the string
alphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
regName = [alphabet{aut.q{iModeToPatch}},'_',alphabet{aut.q{trans(itrans,2)}}];
decompRegName = [];
newStr = [',    {"name":"',regName,'","color":[255,0,0],"holeList":[],"height":0,"points":[['];
for idx = 1:size(vertsToWrite,1)-1
    newStr = [newStr,num2str(vertsToWrite(idx,1)),',',num2str(vertsToWrite(idx,2)),'],['];
end
newStr = [newStr,num2str(vertsToWrite(size(vertsToWrite,1),1)),',',num2str(vertsToWrite(size(vertsToWrite,1),2))];
newStr = [newStr,']],"position":[0.0,0.0],"type":"poly","size":[0.0,0.0]}'];

% write to the file
fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_new.regions','r');
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

fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_new.regions', 'w');
for idx = 1:length(A)
    if A{idx+1} == -1
        fprintf(fid,'%s', A{idx});
        break
    else
        if ~isempty(A{idx})
            if strcmp(A{idx},']')
                fprintf(fid,'%s\n', newStr);
            end
        end
        fprintf(fid,'%s\n', A{idx});
    end
end
fclose(fid);

% write to the file
fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_decomposed_new.regions','r');
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

fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_decomposed_new.regions', 'w');
for idx = 1:length(A)
    if A{idx+1} == -1
        fprintf(fid,'%s', A{idx});
        break
    else
        if ~isempty(A{idx})
            if strcmp(A{idx},']')
                fprintf(fid,'%s\n', newStr);
            end
        end
        fprintf(fid,'%s\n', A{idx});
    end
end
fclose(fid);
