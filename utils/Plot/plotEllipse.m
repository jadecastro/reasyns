function plotEllipse(ell,color)

if dimension(ell(1)) < 3
    for i = 1:length(ell)
        [c,Q] = double(ell(i));
        A = inv(Q);
        Ellipse_plot1(A,c,color);
    end
else
    for i = 1:length(ell)
        [c,Q] = double(ell(i));
        A = inv(Q);
        if isempty(color)
            Ellipse_plot1(A,c);
        else
            plot(ell(i),color);
        end
    end
end