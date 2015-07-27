
global u_sav x_sav x0_sav

reasynsPath

load reasyns_funneldata_22-Mar-2015

pose = ppval(ac_trans{1}.x0.pp,0);

id_region_1 = 1;
id_region_2 = 7;

t_base = 0;
acLastData = {0, 0, 0, 0};

v = 0;
w = 0.001;

ME.message = 'good to go.';
errorMsg = ME.message;

u_sav = [];
x_sav = [];
x0_sav = [];

clk = fix(clock);  % to name the file to save data to
