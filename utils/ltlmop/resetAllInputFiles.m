% Resets the spec, region, structuredslugs, and config files appended with "_new" and "_sav" to their respective pristine states. 

if isunix
    eval(['! cp ',fileName,'.regions ',fileName,'_sav.regions'])
    eval(['! cp ',fileName,'.structuredslugs ',fileName,'_sav.structuredslugs'])
    eval(['! cp ',fileName,'_decomposed.regions ',fileName,'_decomposed_sav.regions'])
    eval(['! cp ',fileName,'.regions ',fileName,'_sav.regions'])
    eval(['! cp ',fileName,'.structuredslugs ',fileName,'_sav.structuredslugs'])
    eval(['! cp ',fileName,'_decomposed.regions ',fileName,'_decomposed_sav.regions'])
    eval(['! cp ',filePath,'/',configAndProblemDomainName,'.config ',filePath,'/',configAndProblemDomainName,'_sav.config '])
elseif ispc
    eval(['! copy ',fileName,'.regions ',fileName,'_sav.regions'])
    eval(['! copy ',fileName,'.structuredslugs ',fileName,'_sav.structuredslugs'])
    eval(['! copy ',fileName,'_decomposed.regions ',fileName,'_decomposed_sav.regions'])
    eval(['! copy ',fileName,'.regions ',fileName,'_sav.regions'])
    eval(['! copy ',fileName,'.structuredslugs ',fileName,'_sav.structuredslugs'])
    eval(['! copy ',fileName,'_decomposed.regions ',fileName,'_decomposed_sav.regions'])
    eval(['! copy ',filePath,'\',configAndProblemDomainName,'.config ',filePath,'/',configAndProblemDomainName,'_sav.config '])
end