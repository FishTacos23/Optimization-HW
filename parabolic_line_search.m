function [a_star, n_objectives] = parabolic_line_search(plot_bool, obj, x, f, s, alpha, mag)
    

    f_a(1) = f(:,end);

    % FIRST STEP
    a(2) = alpha;
    x_n = x(:,end) + a(2)*s(:,end)/mag;
    f_a(2) = obj(x_n);
    
    n_objectives = 1;

    if plot_bool
%         plot(x_n(1), x_n(2), 'r.', 'LineWidth', 5);
    end

    % Iterator Variables
    i = 3;
    not_passed_min = true;
    not_done = true;

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
        x_n = x(:,end) + (a(i)*s(:,end)/mag);
        f_a(i) = obj(x_n);
        
        n_objectives = n_objectives + 1;

        % Iteration funcs
        if plot_bool
%             plot(x_n(1), x_n(2), 'r.', 'LineWidth', 5);
        end

        i = i + 1;

    end

    [val, index2] = min(f_a);
    a_2 = a(index2);
    f_2 = f_a(index2);

    a_less = a;
    a_less(a_less<=a_2)=Inf;
    [val, index3] = min(a_less);
    a_3 = a(index3);
    f_3 = f_a(index3);

    a_more = a;
    a_more(a_more>=a_2)=0;
    [val, index1] = max(a_more);
    a_1 = a(index1);
    f_1 = f_a(index1);

    a_star = (f_1*(a_2^2-a_3^2)+f_2*(a_3^2-a_1^2)+f_3*(a_1^2-a_2^2))/(2*(f_1*(a_2-a_3)+f_2*(a_3-a_1)+f_3*(a_1-a_2)));
    
    
end