
Nskip = 20;
Vol = 2;

sysTrans = sys(2);

acOld = ac_inward{4};
ac = acOld;

for i = 1:5
    ell = projection(ac.ellipsoid,[1 0 0;0 1 0]');
    
    t1 = ac.x0.getTimeVec();
    t2 = t1(volume(ell) < Vol);
    % t2 = downsampleUniformly(t1,Nskip);
    x2 = Traject(t2,double(ac.x0,t2));
    u2 = Traject(t2,double(ac.u0,t2));
    K2 = Traject(t2,double(ac.K,t2));
    P2 = Traject(t2,double(ac.P,t2));
    rho2 = Traject(t2,double(ac.rho,t2));
    % Vquad2 = ac.V(1:Nskip:end);
    Vquad2 = ac.V(volume(ell) < Vol);
    
    % t2 = t1;
    % rhoTmp = double(ac.rho,t2);
    % Ptmp = double(ac.P,t2);
    % VquadTmp = ac.V;
    % iSav = [];
    % for i = 1:length(ell)
    %     if i > 1 && volume(ell(i)) > Vol
    %         iSav = [iSav i];
    %         Ptmp(:,:,i) = Ptmp(:,:,i-1);
    %         VquadTmp(i) = VquadTmp(i-1);
    %         rhoTmp(i) = rhoTmp(i-1);
    %     end
    % end
    % x2 = Traject(t2,double(ac.x0,t2));
    % u2 = Traject(t2,double(ac.u0,t2));
    % K2 = Traject(t2,double(ac.K,t2));
    % rho2 = Traject(t2,rhoTmp);
    % P2 = Traject(t2,Ptmp);
    % Vquad2 = VquadTmp;
    
    acSmoothed = quadraticAC(x2,u2,K2,P2,rho2,Vquad2,sysTrans);
    ac = acSmoothed;
end

% double check
ell2 = projection(acSmoothed.ellipsoid,[1 0 0;0 1 0]');
for i = 1:length(ell2)
    if volume(ell2(i)) > Vol
        i
        disp('it''s baaack!!')
    end
end

% figure(98), clf
% hold on
% ell1 = projection(acSmoothed.ellipsoid,[1 0 0;0 1 0]');
% for i=iSav
%     plot(ell(i),'b')
%     plot(ell1(i),'g')
% end
% return

figure(99), clf
hold on
plot(acOld,sysTrans,99,[],[0 0 1])
plot(acSmoothed,sysTrans,99,[],[0 1 0])
