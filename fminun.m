 function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, algoflag, alpha, num_plot)
             
        %initialize variables
        f(:,1) = obj(x0);          %objective function
        grad(:,1) = gradobj(x0);   %gradient
        x(:,1) = x0;               %design
        s(:,1) = -grad(:,1);
        
        % loop parameters
        iter = 1;          % iteration number
        max_iter = 100;    % maximum number of iterations
        plot_bool = true;
                
        while (iter<max_iter && all((abs(grad(end,:))) > stoptol))
            
            % PLOTTING
            if iter < num_plot
                plot(x(1, end),x(2, end),'cx', 'LineWidth', 10);
            else
                plot_bool = false;
            end
            
            % FIND NEXT STEP
            if (algoflag == 1)     % steepest descent
                [x, f, grad] = steepest_decent(f(end, :), grad(end, :), x(end, :), alpha, obj, gradobj, plot_bool);        
            elseif (algoflag == 2) % conjugate gradient
                [x, f, grad, s] = conjugant_gradient(f, grad, x, alpha, obj, gradobj, plot_bool, s);
            else                   % BFGS quasi newtown
                if (iter==1)
                    N = eye(length(x0));
                end
                [x, f, grad, N] = BFGS(f, grad, x, alpha, obj, gradobj, plot_bool, N);  
            end
            
            iter = iter + 1;        
                        
        end
        
        xopt = x(:,end);
        fopt = f(:,end);
        exitflag = 0;
        
        plot(x(1, end),x(2 ,end),'w*', 'LineWidth', 10);
        
        iter
        
 end