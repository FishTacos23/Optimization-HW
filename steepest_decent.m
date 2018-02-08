 function [x, f, grad, s, a_star, n_objectives] = steepest_decent(f, grad, x, alpha, obj, gradobj, plot_bool, s)

    s(:,end+1) = -grad(:,end);   % search direction
    mag = sqrt(grad(:,end)'*grad(:,end)); % magnitude

    % Set values to calc a_star
    [a_star, n_objectives] = parabolic_line_search(plot_bool, obj, x, f, s, alpha, mag);
    
    x(:,end+1) = x(:,end) + a_star*s(:,end)/mag;
    f(:,end+1) = obj(x(:,end));
    grad(:,end+1) = gradobj(x(:,end));
    
    n_objectives = n_objectives + 1;

 end