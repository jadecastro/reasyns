

figure(2), clf
plot(reg)
hold on
% for i = 1:length(ac_trans)
%     if ~isempty(ac_trans{i})
%         plot(ac_trans{i},sys,2)
%     end
% end
plot(x0_sav(:,2),x0_sav(:,3),'m.')
plot(x_sav(:,2),x_sav(:,3),'k.')
axis equal

figure(3), clf
plot(1:length(x0_sav(:,1)),x0_sav(:,2),'m.',1:length(x_sav(:,1)),x_sav(:,2),'k.')
