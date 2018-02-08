function [x, f, grad, N, s, a_star, n_objectives] = BFGS(f, grad, x, alpha, obj, gradobj, plot_bool, N, s)

    s(:,end+1) = -N*grad(:,end);   % search direction
    mag = sqrt(s(:,end)'*s(:,end)); % magnitude

    [a_star, n_objectives] = parabolic_line_search(plot_bool, obj, x, f, s, alpha, mag);

    x(:,end+1) = x(:, end) + a_star*s(:,end)/mag;
    f(:,end+1) = obj(x(:, end));
    grad(:,end+1) = gradobj(x(:, end));
    
    n_objectives = n_objectives + 1;

    Dx = x(:, end)-x(:,end-1);
    gamma = grad(:, end)-grad(:,end-1);

    N = N + (1+gamma'*N*gamma/(Dx'*gamma))*(Dx*Dx'/(Dx'*gamma))-((Dx*gamma'*N+N*gamma*Dx')/(Dx'*gamma));

 end