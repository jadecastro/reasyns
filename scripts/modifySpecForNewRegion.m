

% compose the string
alphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
newRegName = [alphabet{aut.q{iModeToPatch}},'_',alphabet{aut.q{trans(itrans,2)}}];
currRegName = [alphabet{aut.q{iModeToPatch}}];
nextRegName = [alphabet{aut.q{trans(itrans,2)}}];
preStates = horzcat(aut.q{trans(trans(:,2)==iModeToPatch,1)});
prevRegName = [alphabet{preStates}];

% make sure, when entering R2, we're not already in the IC region
envSafety1 = ['rob1_',currRegName,'_rc & rob1_',currRegName,' & rob1_',currRegName,'_rc'' -> !rob1_',newRegName,'_rc'''];

% conditions for entering the IC region
envSafety2 = ['rob1_',currRegName,'_rc & !rob1_',newRegName,'_rc & rob1_',nextRegName,' -> rob1_',currRegName,'_rc'''];
envSafety3 = ['rob1_',currRegName,'_rc & rob1_',newRegName,'_rc & rob1_',nextRegName,' -> (rob1_',newRegName,'_rc'' | rob1_',nextRegName,'_rc'')'];

% we can always avoid the IC region if not activating R3 or not in R2
for iPrevReg = 1:length(prevRegName)
    envSafety4{iPrevReg} = ['!rob1_',nextRegName,' & m_rob1_',prevRegName(iPrevReg),' -> !rob1_',newRegName,'_rc'''];
end
envSafety5 = ['!rob1_',currRegName,'_rc'' -> !rob1_',newRegName,'_rc'''];

% modify the initial conditions
envInit = ['!rob1_',newRegName,'_rc'];
sysInit = ['!rob1_',newRegName];

% added guarantees
for iPrevReg = 1:length(prevRegName)
    sysSafety1{iPrevReg} = ['m_rob1_',prevRegName(iPrevReg),' & rob1_',currRegName,'_rc'' & rob1_',newRegName,'_rc'' -> rob1_',nextRegName,''''];
    sysSafety2{iPrevReg} = ['((m_rob1_',prevRegName(iPrevReg),' | rob1_',prevRegName(iPrevReg),'_rc'') & (rob1_',prevRegName(iPrevReg),'_rc'' | rob1_',currRegName,'_rc'')) <-> m_rob1_',prevRegName(iPrevReg),''''];
end

% write to the file
fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_new.structuredslugs','r');
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

fid = fopen('/home/jon/Dropbox/Repos/LTLMoP/src/examples/box_pushing/box_pushing_new.structuredslugs', 'w');
for idx = 1:length(A)
    if A{idx+1} == -1
        fprintf(fid,'%s', A{idx});
        break
    else
        
        if strcmp(A{idx},'[OUTPUT]') % TODO: not robust to different orders of the file headings
            fprintf(fid,'%s\n\n', ['rob1_',newRegName,'_rc']);
        end
        if strcmp(A{idx},'[ENV_TRANS]') % TODO: not robust to different orders of the file headings
            fprintf(fid,'%s\n\n', ['rob1_',newRegName,'']);
            fprintf(fid,'%s\n\n', ['m_rob1_',prevRegName,'']);
        end
        
        fprintf(fid,'%s\n', A{idx});
        
        if strcmp(A{idx},'[ENV_TRANS]')
            fprintf(fid,'%s\n', envSafety1);
            fprintf(fid,'%s\n', envSafety2);
            fprintf(fid,'%s\n', envSafety3);
            for iPrevReg = 1:length(prevRegName)
                fprintf(fid,'%s\n', envSafety4{iPrevReg});
            end
            fprintf(fid,'%s\n\n', envSafety5);
        end
        if strcmp(A{idx},'[ENV_INIT]')
            fprintf(fid,'%s\n\n', envInit);
        end
        if strcmp(A{idx},'[SYS_TRANS]')
            for iPrevReg = 1:length(prevRegName)
                fprintf(fid,'%s\n\n', sysSafety1{iPrevReg});
                fprintf(fid,'%s\n\n', sysSafety2{iPrevReg});
            end
        end
        if strcmp(A{idx},'[SYS_INIT]')
            fprintf(fid,'%s\n\n', sysInit);
        end
    end
end
fclose(fid);

% TODO: also modify spec and config files
%
%
%
%



