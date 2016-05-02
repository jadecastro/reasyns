
function saveReasynsProblemInstance(filename,sys,aut,ac_inward,ac_trans,options)

ac_inward_py = [];
k = 0;
for i = 1:length(ac_inward)
    if ~isempty(ac_inward{i})
        for j = 1:length(ac_inward{i})
            k = k + 1;
            ac_inward_py{k} = storeQuadraticAtomicControllerAsStruct(ac_inward{i}(j));
        end
    end
end
k = 0;
for i = 1:length(ac_trans)
    if ~isempty(ac_trans{i})
        for j = 1:length(ac_trans{i})
            k = k + 1;
            ac_trans_py{k} = storeQuadraticAtomicControllerAsStruct(ac_trans{i}(j));
        end
    end
end

eval(['save ',filename,' sys aut ac_inward_py ac_trans_py'])

end

function ac_py = storeQuadraticAtomicControllerAsStruct(ac)

    t = getTimeVec(ac.x0);
    ac_py.t = t;
    ac_py.x0 = double(ac.x0,t);
    ac_py.u0 = double(ac.u0,t);
    ac_py.K = double(ac.K,t);
    ac_py.Einv = ac.Einv;
    ac_py.rho = double(ac.rho,t);
    ac_py.P = double(ac.P,t);
    ac_py.pre = ac.pre;
    ac_py.post = ac.post;

end
