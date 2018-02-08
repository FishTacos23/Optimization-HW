function [x, f, grad, s, a_star, n_objectives] = conjugant_gradient(f, grad, x, alpha, obj, gradobj, plot_bool, s)
    
    mag = sqrt(s(:,end)'*s(:,end)); % magnitude

    [a_star, n_objectives] = parabolic_line_search(plot_bool, obj, x, f, s, alpha, mag);

    x(:,end+1) = x(:,end) + a_star*s(:,end)/mag;
    f(:,end+1) = obj(x(:,end));
    grad(:,end+1) = gradobj(x(:,end));
    
    n_objectives = n_objectives + 1;
    
    Beta = grad(:,end)'*grad(:,end)/(grad(:,end-1)'*grad(:,end-1));
    s(:,end+1) = -grad(:,end) + Beta*s(:,end);

 end