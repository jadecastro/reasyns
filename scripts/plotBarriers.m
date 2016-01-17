
function plotBarriers(ac_tmp,bc_tmp)

N = length(bc_tmp);
[uRed,tRed] = downsampleUniformly(ac_tmp.x0,10);
colorArray = colormap('cool');
x = msspoly('x',3);

for ii = 1:N
    color = colorArray(floor(ii*size(colorArray,1)/N),:);
    
%     prog(B{ii})

    % check for sanity- check equilibrium point.  TODO: check more points in
    % the IC set too

    x0 = double(ac_tmp.x0,tRed(ii));
    
    % 0-level set of B
    [X,Y] = meshgrid(linspace(-3,2,100),linspace(-2,2,100));
    Th = x0(3)*ones(size(X));
    gs_B = msubs(bc_tmp{ii},x,[X(:),Y(:),Th(:)]');
    [~,H] = contour(X,Y,reshape(double(gs_B),100,100),[0 0],'LineColor',color,'LineWidth',3); 
    %set(H,'LineColor',color)
    hold on
    
%     % testing...
%     [ellq,ellQ] = double(tmpArray(idx00));
%     g_test = 0.01 - (x - ellq)'*inv(ellQ)*(x - ellq);
%     Th = ellq(3)*ones(size(X));
%     gs_test = msubs(sol.eval(g_test),x,[X(:),Y(:),Th(:)]');
%     [~,H] = contour(X,Y,reshape(double(gs_test),100,100),[0 0],'LineColor',color,'LineWidth',3); 

end