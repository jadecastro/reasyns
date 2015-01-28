
for i = 1:5, 
    Ktk{1}{i} = funnel{i}.K; 
    utk{1}{i} = funnel{i}.u; 
    xtk{1}{i} = funnel{i}.x; 
    ttk{1}{i} = funnel{i}.t; 
    Ptk{1}{i} = funnel{i}.P; 
end
for i = 1:4, 
    for j = 1:5, 
        Kik{j}{i} = funnelIn{j,i}.K; 
        uik{j}{i} = funnelIn{j,i}.u; 
        xik{j}{i} = funnelIn{j,i}.x; 
        tik{j}{i} = funnelIn{j,i}.t; 
        Pik{j}{i} = funnelIn{j,i}.P; 
    end, 
end
for i = 1:5, 
    Kjk{1}{i} = funnelJoin{i}.K; 
    ujk{1}{i} = funnelJoin{i}.u; 
    xjk{1}{i} = funnelJoin{i}.x; 
    tjk{1}{i} = funnelJoin{i}.t; 
    Pjk{1}{i} = funnelJoin{i}.P; 
end

% TODO: reactive join

save('~/testwithstruct.mat','ttk','xtk','utk','Ktk','Ptk','tik','xik','uik','Kik','Pik','tjk','xjk','ujk','Kjk','Pjk','aut')