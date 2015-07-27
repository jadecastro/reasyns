function writeToFlowstarModelFile(infname, outfname, xi_box, xf)
% Tranfer data from a temlpate data file to another, replacing any
% 'placeholders' with custom data.

fid = fopen(infname,'r');
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

% construct the needed strings
initBoxString = sprintf('x in [%s,%s] y in [%s,%s] theta in [%s,%s]\n',num2str(xi_box(1,1)),num2str(xi_box(1,2)),num2str(xi_box(2,1)),num2str(xi_box(2,2)),num2str(xi_box(3,1)),num2str(xi_box(3,2)));
invString = sprintf('inv { (x - (%s))^2 + (y - (%s))^2 - 0.25^2 >= 0}\n',num2str(xf(1)),num2str(xf(2)));
thetaDotString = sprintf('theta'' = 1/0.1*(-(%s - x)*sin(theta) + (%s - y)*cos(theta))\n',num2str(xf(1)),num2str(xf(2)));

fid = fopen(outfname, 'w');
for idx = 1:length(A)
    placeholderFlag = false;
    if A{idx+1} == -1
        fprintf(fid,'%s', A{idx});
        break
    else
        if strfind(A{idx},'init_placeholder') 
            placeholderFlag = true;
            fprintf(fid,initBoxString);
        end
        if strfind(A{idx},'inv_placeholder') 
            placeholderFlag = true;
            fprintf(fid,invString);
        end
        if strfind(A{idx},'thetadot_placeholder')
            placeholderFlag = true;
            fprintf(fid,thetaDotString);
        end
        if ~placeholderFlag
            fprintf(fid,'%s\n', A{idx});
        end
    end
end
fclose(fid);