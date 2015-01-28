
load testCase10

ellFunnel0 = ellFunnel;
trajFunnel0 = trajFunnel;
uFunnel0 = uFunnel;
Pfunnel0 = Pfunnel;
Kfunnel0 = Kfunnel;
clear ellFunnel trajFunnel uFunnel Pfunnel Kfunnel
for j = 1:size(ellFunnel0,2)
    for i = 1:size(ellFunnel0,1)
        if ~isempty(ellFunnel0{i,j,1})
%             elim(i) = max(maxeig(ellFunnel0{i,j,1})) < 0.05;
            elim(i,j) = rhoMinArray{i,j,1} < 0.5;
            if i > 1 && (sum(any(Isect{j,1}(1:i,:),1)) - sum(any(Isect{j,1}(1:i-1,:),1)) == 0)
                elim(i,j) = 1;
            end
        end
    end
end
for j = 1:size(ellFunnel0,2)    
    count = 0;
    for i = 1:size(ellFunnel0,1)
        if ~elim(i,j),
            count = count+1;
            ellFunnel{count,j,1} = ellFunnel0{i,j,1};
            trajFunnel{count,j,1} = trajFunnel0{i,j,1};
            uFunnel{count,j,1} = uFunnel0{i,j,1};
            Pfunnel{count,j,1} = Pfunnel0{i,j,1};
            Kfunnel{count,j,1} = Kfunnel0{i,j,1};
        end
    end
end

% ellFunnel_tmp = ellFunnel;
% trajFunnel_tmp = trajFunnel;


%%

load testCase40

ellFunnelC0 = ellFunnelC;
trajFunnelC0 = trajFunnelC;
uFunnelC0 = uFunnelC;
PfunnelC0 = PfunnelC;
KfunnelC0 = KfunnelC;
clear ellFunnelC trajFunnelC uFunnelC PfunnelC KfunnelC
elimC = zeros(size(ellFunnelC0,1),size(ellFunnelC0,2));
for j = 1:size(ellFunnelC0,2)
    for i = 1:size(ellFunnelC0,1)
        if ~isempty(ellFunnelC0{i,j,2})
%             elim(i,j) = max(maxeig(ellFunnelC0{i,j,1})) < 0.05;
            elimC(i,j) = rhoMinArrayC{i,j,2} < 0.05;
        end
    end
end
for j = 1:size(ellFunnelC0,2)    
    count = 0;
    for i = 1:size(ellFunnelC0,1)
        if ~elimC(i,j),
            count = count+1;
            ellFunnelC{count,j,2} = ellFunnelC0{i,j,2};
            trajFunnelC{count,j,2} = trajFunnelC0{i,j,2};
            uFunnelC{count,j,2} = uFunnelC0{i,j,2};
            PfunnelC{count,j,2} = PfunnelC0{i,j,2};
            KfunnelC{count,j,2} = KfunnelC0{i,j,2};
        end
    end
end

%%
figure(7)
clf
hold on

options.fill = 1;
options.color = [0.9 0.3 0.3];

for j = 1:size(ellFunnel,2)
    for i = 1:size(ellFunnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plot(ellFunnel{i,j}(1)); end
        if ~isempty(ellFunnel{i,j,2}), plotEllipse(projection(ellFunnel{i,j,2}(1),[1 0;0 0;0 1]),options.color); end
    end
end

options.color = [0.3 0.9 0.3];

for j = 1:size(ellFunnelC,2)
    for i = 1:size(ellFunnelC,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(ellFunnelC{i,j}(1),'g'); end
        if ~isempty(ellFunnelC{i,j,2}), plotEllipse(projection(ellFunnelC{i,j,2}(1),[1 0;0 0;0 1]),options.color); end
    end
end


H(1) = hyperplane([1 0 0]',0.1);
H(2) = hyperplane([1 0 0]',0.5);
H(3) = hyperplane([1 0 0]',1);
H(4) = hyperplane([1 0 0]',1.75);
H(5) = hyperplane([1 0 0]',2);
H(6) = hyperplane([1 0 0]',2.5);

options.fill = 1;
options.color = [0.9 0.3 0.3];

for ii = 1:length(H)
    figure(ii+7)
    clf
    hold on
    for j = 1:size(ellFunnel,2)
        for i = 1:size(ellFunnel,1)
            %         if ~isempty(ellFunnel{i,j,1}), plot(ellFunnel{i,j}(1)); end
            if ~isempty(ellFunnel{i,j,2})
                tmp = hpintersection(ellFunnel{i,j,2}(1),H(ii));
                if ~isempty(tmp), plotEllipse1(tmp,options.color); end
            end
        end
    end
end

options.color = [0.3 0.9 0.3];

for ii = 1:length(H)
    figure(ii+7)
    for j = 1:size(ellFunnelC,2)
        for i = 1:size(ellFunnelC,1)
            %         if ~isempty(ellFunnelC{i,j,1}), plot(ellFunnelC{i,j}(1),'g'); end
            if ~isempty(ellFunnelC{i,j,2})
                tmp = hpintersection(ellFunnelC{i,j,2}(1),H(ii));
                if ~isempty(tmp), plotEllipse1(tmp,options.color); end
            end
        end
    end
end

%%

figure(3)
clf
hold on
axis equal

options.color = [0.9 0.9 0.9];
plot(pReg{2},options)
plot(pReg{1},options)
plot(pReg{3},options)
figure(4)
clf
hold on
axis equal
plot(pReg{2},options)
plot(pReg{1},options)
plot(pReg{3},options)

options.fill = 1;
options.color = [0.9 0.3 0.3];
for j = 1:size(ellFunnel,2)
    for i = 1:size(ellFunnel,1)
        if ~isempty(ellFunnel{i,j,2}), plotEllipse(projection(ellFunnel{i,j,2},[1 0;0 1;0 0]),options.color); end
    end
end
for j = 1:size(ellFunnel,2)
    for i = 1:size(ellFunnel,1)
        if ~isempty(ellFunnel{i,j,2}), plot(trajFunnel{i,j,2}(:,1),trajFunnel{i,j,2}(:,2),'k','LineWidth',2); end
    end
end

options.fill = 1;
options.color = [0.3 0.9 0.3];
for j = 1:size(ellFunnelC,2)
    for i = 1:size(ellFunnelC,1)
        if ~isempty(ellFunnelC{i,j,2}), plotEllipse(projection(ellFunnelC{i,j,2},[1 0;0 1;0 0]),options.color); end
    end
end
for j = 1:size(ellFunnelC,2)
    for i = 1:size(ellFunnelC,1)
        if ~isempty(ellFunnelC{i,j,2}), plot(trajFunnelC{i,j,2}(:,1),trajFunnelC{i,j,2}(:,2),'k','LineWidth',2); end
    end
end