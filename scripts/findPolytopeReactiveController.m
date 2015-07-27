
% load a funnel
load('reasyns_data_14-Feb-2015_185727.mat')


% make a region for the initial set
clear reg0
reg0 = region([X{1}(1:2:end-1),Y{1}(1:2:end-1)]);

figure(1)
clf
plot(reg)
axis equal
hold on 

plot(ac_trans{3},sys,1,0)

plot(reg0)

%%
% compute alignmentAngle
%TODO: compute hyperplane index from intersection of the transition reach
%set with a particular facet
%TODO: orient alignmentAngle depending on the intersection of the init set with the
%reach set
%TODO: decide between alignmentAngle or alignmentAngle+pi based on containment in the region
idx = 4;        % hyperplane index
polyIdx = 2;    % index of the polytope
[Hp,Kp] = double(reg(1).p(polyIdx));
alignmentAngle = atan2(-Hp(idx,1)/Hp(idx,2),1);
alignmentAngle = alignmentAngle + pi/2;  % convert from tangent to hyperplane normal
alignmentAngle = alignmentAngle + pi;  % TODO^ 

for i = 1:length(X)
    i
    rotVert{i} = [X{i}(1:end-1) Y{i}(1:end-1)]*[cos(alignmentAngle), -sin(alignmentAngle); sin(alignmentAngle), cos(alignmentAngle)]';
end
reg0 = region(rotVert{57});
plot(reg0)

%%
% find points that are on the line
[Hp0,Kp0] = double(reg0.p);
for origIdx = 1:size(Hp0,1)
    for lineIdx = 1:size(Hp,1)
        if norm(Hp(lineIdx,:)-Hp0(origIdx,:)) < 1e-8
            break
        end
    end
    if norm(Hp(lineIdx,:)-Hp0(origIdx,:)) < 1e-8
        break
    end
end

% what other lines in the shape are parallel?
for origIdxOpp = 1:size(Hp,1)
    if origIdxOpp~=origIdx && norm(abs(Hp(lineIdx,:))-abs(Hp0(origIdxOpp,:))) < 1e-8
        break
    end
end

[~,Idx] = unique(reg0.v(:,1));
X0 = reg0.v(sort(Idx),1);  Y0 = reg0.v(sort(Idx),2);

kp_delta = Kp(lineIdx) - Kp0(origIdx);
Xnew = X0;% + kp_delta/Hp0(origIdx,1);
Ynew = Y0 + kp_delta/Hp0(origIdx,2);

reg1 = region([Xnew Ynew]);

hold on
plot(reg1)

%%
tangIdx = 2;

%TODO: automate the finding of the multiplier based on a greedy search
multiplier = -0.15;
Xnew = Xnew + multiplier*kp_delta/Hp0(tangIdx,2);
Ynew = Ynew + multiplier*kp_delta/Hp0(tangIdx,1);

% collect the transformation offsets here
Xoff = multiplier*kp_delta/Hp0(tangIdx,2);
Yoff = kp_delta/Hp0(origIdx,2) + multiplier*kp_delta/Hp0(tangIdx,1);

reg2 = region([Xnew Ynew]);

hold on
plot(reg2)
for i = 1:length(X)
    if length(X) < 1000 || ~rem(i,10)
        plot([rotVert{i}(:,1); rotVert{i}(1,1)]+Xoff,[rotVert{i}(:,2); rotVert{i}(1,2)]+Yoff,'Color',[0.5 0.5 0.5])
    end
end

plot(reg2)

%%
reg00 = region(rotVert{1});
[Hp00,Kp00] = double(reg00.p);

[~,Idx] = unique(reg00.v(:,1));
X00 = reg00.v(sort(Idx),1);  Y00 = reg00.v(sort(Idx),2);

Xnew = X00 + multiplier*kp_delta/Hp0(tangIdx,2);
Ynew = Y00 + kp_delta/Hp0(origIdx,2) + multiplier*kp_delta/Hp0(tangIdx,1);

reg00 = region([Xnew Ynew]);

plot(reg00)

