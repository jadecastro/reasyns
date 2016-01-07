
global u_sav x_sav x0_sav

reasynsPath

% "Complex Map" scenario:
%load reasyns_funneldata_22-Mar-2015

% pose = ppval(ac_trans{1}.x0.pp,0);
% 
% id_region_1 = 1;
% id_region_2 = 7;


% "Box Pushing" scenario:
%load box_pushing_funnels_14-Aug-2015
%load box_pushing_funnels_13-Sep-2015b_oneInwardFunnel

%load box_pushing_lab_funnels_19-Aug-2015
%load box_pushing_lab_funnels_27-Sep-2015_oneInwardFunnel
%load box_pushing_lab_funnels_29-Sep-2015_oneInwardFunnel
load box_pushing_lab_funnels_29-Sep-2015

clear ME

pose = ppval(ac_trans{1}.x0.pp,0);

id_region_1 = 5;
id_region_2 = 1;


t_base = 0;
acLastData = {1, 0, 5, 1, [], []};

v = 0;
w = 0.001;

ME.message = 'good to go.';
errorMsg = ME.message;

u_sav = [];
x_sav = [];
x0_sav = [];

clk = fix(clock);  % to name the file to save data to

global u_sav x_sav x0_sav
