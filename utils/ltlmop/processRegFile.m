
function [reg,regDefl,regBnd] = processRegFile(fpath,fname,options,varargin)
% Load the regions from a .regions file in the required LTLMoP format.

debug = true;

if nargin > 3
    pix2m = varargin{1};
else
    pix2m = 1;
end

disp('Loading the region file...')

% First, get the region and config file names
specData = readFile([fpath,'/',fname,'.spec']);

regKey   = 'RegionFile';
cfgKey   = 'CurrentConfigName';

regfname = [];  configfname = [];
for idx = 1:length(specData)
    regIndex = strfind(specData{idx}, regKey);
    cfgIndex = strfind(specData{idx}, cfgKey);
    if ~isempty(regIndex)
        i1 = idx;
        while isempty(regfname) || regfname(1) == '#'
            i1 = i1+1;
            regfname = specData{i1};
        end
    end
    if ~isempty(cfgIndex)
        i2 = idx;
        while isempty(configfname) || configfname(1) == '#'
            i2 = i2+1;
            configfname = specData{i2};
        end
    end
end

% Retrieve and process the region file data structure
fid = fopen([fpath,'/',regfname]);

% We need to extract the json part of the data file first
regFileString = [];
flag = 0;
while true
    line = fgetl(fid);
    
    if ~ischar(line), regFileString = [regFileString ' ]']; break, end
    
    if strfind(line,'[')
        flag = flag + 1;
    end
    if strfind(line,']')
        flag = flag - 1;
    end
    if flag
        regFileString = [regFileString line];
    end
end
fclose(fid);

% Parse the json part into a data structure
regData = loadjson(regFileString);

% Now, retrieve the calibration matrix
cfgData = readFile([fpath,'/configs/',configfname,'.config']);

calKey   = 'CalibrationMatrix';
endKey   = ')';

% NB: this implementation requires that each row of the calibration matrix
% contains exactly 3 entries in the config file.
calibMatrix = [];
for idx = 1:length(cfgData)
    calIndex = strfind(cfgData{idx}, calKey);
    if ~isempty(calIndex)
        i = idx;
        while isempty(calibMatrix)
            i = i+1;
            if isempty(strfind(cfgData{i},'#') == 1) % skip this line if commented
                tmpStringArray = regexp(cfgData{i},'(?<=\s)[-\d.]+(?=\D)','match');
                if length(tmpStringArray) ~= 3, error('Calibration matrix must have exactly 3 columns per row.'); end
                for idx = 1:3, tmp(idx) = str2num(tmpStringArray{idx}); end
                calibMatrix = tmp;
                
                % handle multi-line calibration matrix entries
                if ~isempty(calibMatrix)
                    endIndex = strfind(cfgData{i}, endKey);
                    j = i;
                    while isempty(endIndex)
                        j = j+1;
                        endIndex = strfind(cfgData{j}, endKey);
                        
                        tmpStringArray = regexp(cfgData{j},'(?<=\s)[-\d.]+(?=\D)','match');
                        if length(tmpStringArray) ~= 3, error('Calibration matrix must have exactly 3 columns per row.'); end
                        for idx = 1:3, tmp(idx) = str2num(tmpStringArray{idx}); end
                        calibMatrix = [calibMatrix; tmp];
                    end
                end
                if size(calibMatrix,1) ~= 3, error('Calibration matrix must have exactly 3 rows.'); end
            end
        end
    end
end
               
% Construct the regions
idx = 1;
for i = 1:length(regData)
    if regData{i}.type == 'rect'
        vert = [];
        vert(1,:) = regData{i}.position;
        vert(2,:) = regData{i}.position + [regData{i}.size(1) 0];
        vert(3,:) = regData{i}.position + [regData{i}.size(1) regData{i}.size(2)];
        vert(4,:) = regData{i}.position + [0 regData{i}.size(2)];
        
    elseif regData{i}.type == 'poly'
        numPoints = size(regData{i}.points,1);
        vert = repmat(regData{i}.position,numPoints,1) + regData{i}.points;
        
    else
        error('Unrecognized region type.')
    end
        
    if strmatch(regData{i}.name,'boundary')
        regBnd = Region(regData{i}.name, vert, calibMatrix);  % TODO: deprecate regBnd
    else
        reg(idx)     = Region(regData{i}.name, vert, calibMatrix);
        regDefl(idx) = Region(regData{i}.name, vert, calibMatrix, -options.deflationAmount);  
        idx = idx+1;
    end
end

if debug
    figure, plot(reg)
end

% Perform a sanity check on the loaded regions: 
% check that we have unique region names and that each region has a nonzero intersection with other regions.
disp('Checking the loaded regions...')
for k = 1:length(reg)
    regInfl(k) = Region(reg(k).name,reg(k).v,1, 0.01);
end

for i = 1:length(reg)
    allNames{i} = reg(i).name;
    
    % pairwise check each region
    emptyIntersects = [];
    for j = i:length(reg)
        emptyIntersects(j) = isempty(intersect(regInfl(i),regInfl(j)));
    end
    if sum(~emptyIntersects) < 1, error('Region file must contain overlapping regions!'); end
end

uniqueRegionNames = unique(allNames(:));

if length(uniqueRegionNames) ~= length(reg), error('Region file must contain unique region identifiers!'); end


function [A] = readFile(fname)
fid = fopen(fname,'r');

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
