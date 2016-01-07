

figure(2), clf
plot(reg(1:5))
hold on
% for i = 1:length(ac_trans)
%     if ~isempty(ac_trans{i})
%         plot(ac_trans{i},sys,2)
%     end
% end

Nstart = 250;
Nend = length(x_sav);
% Nend = 255;
%Nend = 275;
% Nend = 350;
% Nend = 500;

plot(reg(6))
H = get(get(gco,'Children'),'Children');
%set(H(1),'FaceColor',[1,0.5,0.5])
HA = get(gco,'Children');
% set(HA,'XTick',[],'XLabel',[],'YTick',[],'YLabel',[])
set(HA,'XTick',[],'YTick',[])
xlabel(''); ylabel('');

% plot funnel
plot(ac_trans{1},sys(1),2,[],[0.6 0.6 0.8])
% plot(ac_trans{4},sys,2,[],[0.8 0.6 0.6])
% plot(acPost,sys,2,[],[0.6 0.8 0.6])

%plot(x0_sav(Nstart:Nend,2),x0_sav(Nstart:Nend,3),'LineWidth',2)
plot(x0_sav(Nstart:Nend,2),x0_sav(Nstart:Nend,3),'m.','LineWidth',2)
%plot(x_sav(Nstart:Nend,2),x_sav(Nstart:Nend,3),'k','LineWidth',4,'Color',[0.5,0.5,0.5])
plot(x_sav(Nstart:Nend,2),x_sav(Nstart:Nend,3),'Color',[0.1,0.1,0.1],'LineWidth',4)
axis equal

figure(3), clf
plot(1:length(x0_sav(:,1)),x0_sav(:,2),'m.',1:length(x_sav(:,1)),x_sav(:,2),'k.')

figure(2)
for i = Nstart:55:Nend
    %drawCar(x_sav(i,2:4)')
end
    
