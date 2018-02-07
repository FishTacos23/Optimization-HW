function [s, a_star, mag] = conjugant_gradient(f, grad, x, alpha, obj, plot_bool)

    s = -grad;   % search direction
    mag = sqrt(grad'*grad); % magnitude

    f_a(1) = f;

    % FIRST STEP
    a(2) = alpha;
    x_n = x + a(2)*s/mag;
    f_a(2) = obj(x_n);

    if plot_bool
        plot(x_n(1), x_n(2), '--r*');
    end

    % Iterator Variables
    i = 3;
    not_passed_min = true;
    not_done = true;

    while not_done

        % Decide to step forward or backward
        if not_passed_min == false
            a(i) = a(i-1)/2;
        elseif f_a(i-1) < f_a(i-2)
            a(i) = a(i-1)*2;
        else
            a(i) = a(i-1)*3/4;
            not_passed_min = false;
        end

        [val, min_ind] = min(f_a);

        % Determine if ready for quad
        if not_passed_min == false && min_ind ~= 1
            not_done = false;
        end

        % set point
        x_n = x + a(i)*s/mag;
        f_a(i) = obj(x_n);

        % Iteration funcs
        if plot_bool
            plot(x_n(1), x_n(2), 'r*');
        end

        i = i + 1;

    end

    % Set values to calc a_star       
    a_star = parabolic_line_search(a, f_a);

 end