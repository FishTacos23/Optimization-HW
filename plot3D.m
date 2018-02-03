function [] = plot3D(x1_l,x1_u,x2_l,x2_u,x3)

    % design variables at mesh points
    [x1,x2] = meshgrid(x1_l:(x1_u-x1_l)/500:x1_u,x2_l:(x2_u-x2_l)/500:x2_u);
    
    f = 20+ 3*x1-6*x2+8*x3+6*x1.^2-2*x1.*x2-x1.*x3+x2.^2+0.5*x3.^2;
    
    figure(1)
    contourf(x1,x2,f,-30:10:600);
    title('Quadratic');
    xlabel('x1');
    ylabel('x2');
    hold on;
       
end