function plotEllipse1(ell,color,flg,transparency)

[~,n2] = dimension(ell);

if nargin < 4
    transparency = 1;
end

for i = 1:length(ell)
    
    if n2 == 2
        if strcmp(flg,'waterfall')
            [c,Q] = double(ell(i));
            A = inv(Q(diag(Q)~=0,diag(Q)~=0));
            [X,Y] = Ellipse_plot1(A,c(diag(Q)~=0),color);
            tmp(:,1) = X';
            tmp(:,2) = Y';
            XX(1:length(X),diag(Q)==0) = c(diag(Q)==0);
            count = 0;
            for j = find(diag(Q)~=0)'
                count = count+1;
                XX(:,j) = tmp(:,count);
            end
            H = fill3(XX(:,1),XX(:,2),XX(:,3),color);
            set(H,'LineStyle','none','FaceColor',color,'FaceAlpha',1)

        elseif strcmp(flg,'2d')
            [c,Q] = double(ell(i));
            A = inv(Q(diag(Q)~=0,diag(Q)~=0));
            [X,Y] = Ellipse_plot1(A,c(diag(Q)~=0),color);
            XX(:,1) = X';
            XX(:,2) = Y';
            H = fill(XX(:,1),XX(:,2),color);
            set(H,'LineStyle','none','FaceColor',color,'FaceAlpha',transparency)
                   
        elseif strcmp(flg,'2d_outline')
            [c,Q] = double(ell(i));
            A = inv(Q(diag(Q)~=0,diag(Q)~=0));
            [X,Y] = Ellipse_plot1(A,c(diag(Q)~=0),color);
            XX(:,1) = X';
            XX(:,2) = Y';
            H = fill3(XX(:,1),XX(:,2),-ones(size(XX(:,1))),color);
            set(H,'LineStyle','-','LineWidth',3,'FaceColor','none')
             
        end
        
    elseif n2 == 3
        [c,Q] = double(ell(i));
        A = inv(Q);
        Ellipse_plot1(A,c,color);
    end
end