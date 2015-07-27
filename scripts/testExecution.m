
load reasyns_data_08-Feb-2015_13855
x_sav1 = x_sav;
initExecute
acLastData = {1, 0,1, 2,[],[]};

for i = 54:60
    pose = x_sav1(i,2:4);
    [v,w,errorMsg,acLastData] = executeControllersSingleStep(aut,ac_trans,ac_inward,pose,[],7,1,acLastData)
end