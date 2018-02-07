 function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, algoflag, alpha, num_plot)
             
        %initialize variables
        f(1,:) = obj(x0);          %objective function
        grad(1,:) = gradobj(x0);   %gradient
        x(1,:) = x0;               %design
        
        % loop parameters
        iter = 1;          % iteration number
        max_iter = 100;    % maximum number of iterations
        plot_bool = true;
                
        while (iter<max_iter && all((abs(grad(end,:))) > stoptol))
            
            % PLOTTING
            if iter < num_plot
                plot(x(end, 1),x(end, 2),'-wp');
            else
                plot_bool = false;
            end
            
            % FIND NEXT STEP
            if (algoflag == 1)     % steepest descent
                [x, f, grad] = steepest_decent(f(end, :), grad(end, :), x(end, :), alpha, obj, gradobj, plot_bool);        
            elseif (algoflag == 2) % conjugate gradient
                [x, f, grad] = conjugant_gradient(f(end, :), grad(iter, :), x(iter, :), alpha, obj, plot_bool);
            else                   % BFGS quasi newtown
                if (iter==1)
                    N = eye(length(x0));
                end
                [x, f, grad, N] = BFGS(f, grad, x, alpha, obj, gradobj, plot_bool, N);  
            end
            
            iter = iter + 1;        
                        
        end
        
        xopt = x(end, :);
        fopt = f(end, :);
        exitflag = 0;
        iter
        
 end