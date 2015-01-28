                  for trial1 = 1:maxTrials1
                        [randPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,[],[],'discrep',[],hindx);  % Get initial point
                        initXYTh = randPt;
                        initXY = randPt(1:2);
                        finalPt = getCenter(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,[]);  % Get final point
                        finalPt
                        goalXY = finalPt(1:2);
                        try
                            disp('Building RRT....')
                            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                            [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                            isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                        catch error
                            rethrow(error)
                            isect = true;
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                        for trial2 = 1:maxTrials2
                            disp('Trajectory incompatible with constraints; recomputing...')
                            
                            [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,[],[],'discrep',[],hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                            finalPt
                            goalXY = finalPt(1:2);
                            try
                                disp('Building RRT....')
                                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                            catch error
                                rethrow(error)
                                isect = true;
                            end
                            if ~isect && length(t) > 1 && length(t) < 300, break, end
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                    end