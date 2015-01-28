function plotWS(pReg)

color = {'g','b','y','b','y','b','y','b','y',[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1],[0.1 0.1 0.1]};
Options.shade = 0.5;
for i = 3:5
    figure(i)
    clf
    hold on
    axis equal
    for j = 1:length(pReg)
        Options.color = color{j};
        plot(pReg{j},Options)
    end
end