 function [xopt, fopt, exitflag] = fminun(obj, gradobj, x0, stoptol, ...
     algoflag, alpha, num_plot)
             
        %initialize variables
        f(:,1) = obj(x0);          %objective function
        grad(:,1) = gradobj(x0);   %gradient
        x(:,1) = x0;               %design
        s(:,1) = -grad(:,1);
        alphaPrime = (0);
        nobj = (0);
        
        % loop parameters
        iter = 1;          % iteration number
        max_iter = 500;    % maximum number of iterations
        plot_bool = true;
        
        % SAVE INFORMATION
        saveTable = table;
        mag = sqrt(s(:,end)'*s(:,end));
        dir = s(:,end)'/mag;
        newRow = {x(:,end)', f(:,end)', dir, alphaPrime(:,end)', ...
            nobj(:,end)'};
        saveTable = [saveTable; newRow];
                       
        while (iter<max_iter && all((abs(grad(end,:))) > stoptol))
                         
            % PLOTTING
            if iter < num_plot
                plot(x(1, end),x(2, end),'cx', 'LineWidth', 5);
            else
                plot_bool = false;
            end
            
            % FIND NEXT STEP
            if (algoflag == 1)     % steepest descent
                [x, f, grad, s, a_star, n_step] = steepest_decent(f, ...
                    grad, x, alpha, obj, gradobj, plot_bool, s);        
            elseif (algoflag == 2) % conjugate gradient
                [x, f, grad, s, a_star, n_step] = conjugant_gradient(f,...
                    grad, x, alpha, obj, gradobj, plot_bool, s);
            else                   % BFGS quasi newtown
                if (iter==1)
                    N = eye(length(x0));
                end
                [x, f, grad, N, s, a_star, n_step] = BFGS(f, grad, x,...
                    alpha, obj, gradobj, plot_bool, N, s);  
            end
            
            iter = iter + 1;
            alphaPrime(end+1) = a_star;
            nobj(end+1) = n_step;
            
            mag = sqrt(s(:,end)'*s(:,end));
            dir = s(:,end)'/mag;
            
            newRow = {x(:,end)', f(:,end)', dir, alphaPrime(:,end)',...
                nobj(:,end)'};
            saveTable = [saveTable; newRow];
                        
        end
        
        xopt = x(:,end);
        fopt = f(:,end);
        exitflag = 0;
        
        plot(x(1, end),x(2 ,end),'g*', 'LineWidth', 10);
        
        x1_line = x(1,:)';
        x2_line = x(2,:)';
        plot(x1_line, x2_line, 'k', 'LineWidth', 2);
        
        iter
        
        toSave = table2array(saveTable);
        fout = fopen(sprintf('output%d.csv', algoflag),'w');
        fprintf(fout, '%s, %s, %s, %s, %s, %s, %s, %s, %s\r\n', 'a',...
            'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i');
        fprintf(fout, '%8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f\r\n',...
            toSave');
     	fclose(fout);
        
 end