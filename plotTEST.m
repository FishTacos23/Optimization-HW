function [] = plotTEST(x1_l,x1_u,x2_l,x2_u)

    % design variables at mesh points
    [x1,x2] = meshgrid(x1_l:(x1_u-x1_l)/500:x1_u,x2_l:(x2_u-x2_l)/500:x2_u);
    
    f = x1.^2 - 2*x1.*x2 + 4*x2.^2;
    
    figure(1)
    contourf(x1,x2,f,-1:1:40);
    title('Quadratic');
    xlabel('x1');
    ylabel('x2');
    hold on;
       
end