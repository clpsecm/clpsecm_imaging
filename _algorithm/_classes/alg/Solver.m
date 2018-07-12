classdef Solver < handle
% SOLVER is an abstract framework of algorthms that solves optimization 
% problem using iteration based numerical method.
%
% obj = SOLVER(vars) Construct solver with vars as cell of variables.
% 
% SOLVER is a handle object, and is an abstract class.
%
% SOLVER methods:
%   SET_DISPLAY    - Enable/Disable display and setup frequency of display.
%   SET_MAXITER    - Set maximum iteration number for algorithm.    
%   SOLVE          - Solve problem with implemented algorithm.
%   DISPLAY_RESULT - Show the status of current iteration.
%
% SOLVER abstract methods:
%   ITER           - One iteration of algorithm.
%   SET_OBJVAL     - Calculate the objective value.
%   STOP_CRITERIAS - Evaluate if algorithm satisfies stopping conditions.
%
% SOLVER public fields:
%   OBJVAL  - scalar; objective value.
%   NVARS   - scalar; number of variables.
%   VARS    - cell(nvar,1); variables at current iterates.
%   IITER   - scalar; current iteration number.
%   STATUS  - char; current solver condition.
%

properties
    objval  % scalar; objective value.
    nvars   % scalar; number of variables.
    vars    % cell(nvar,1); variables at current iterates.
    iiter   % scalar; current iteration number.
    status  % char; current solver condition.
end

properties (Hidden)
    objvals_      % vector(maxiter,1) past objective value.
    maxiter       % scalar; maximum number of iterations.
    niter_display % scalar; number of iteration to demo result.
    plot_on       % binary; plot on/off in show_iters().
end

methods (Abstract, Access = protected)
    iter(obj);                      % one iteration of algorithm
    set_objval(obj);                % set objective value
    stopping = stop_criterias(obj); % stopping criteria 
end

methods
    function obj = Solver(vars)
    % SOLVER is a contructor of general framework of iterative solver.
        obj.objval = inf;
        obj.vars  = vars;
        obj.nvars = length(obj.vars);
        obj.iiter = 0;
        obj.status = 'Unsolved';
        obj.objvals_ = [obj.objval; nan(obj.maxiter-1,1)];
        obj.set_maxiter(5000);
        obj.set_display(20);
    end

    function solve(obj)
    % obj.SOLVE Solves the problem using implemented method.
        obj.display_result();
        while true
            obj.iiter = obj.iiter + 1;
            obj.iter();
            obj.objvals_(obj.iiter) = obj.objval;
            if obj.stop(); break; end
            obj.display_iters();
        end
        obj.set_objval();
        obj.display_result();
    end

    function set_maxiter(obj,maxiter); obj.maxiter = maxiter; end
    % obj.SET_MAXITER(maxiter) Set maximum iteration number.

    function set_display(obj, niter_display)
    % obj.SET_DISPLAY(niter_display) Set display setting with frequency
    % niter_display. Input 0 to disable plot.
        obj.plot_on = (niter_display ~= 0);
        obj.niter_display = niter_display;
    end

    function display_result(obj)
    % obj.DISPLAY_RESULT Show the status of current iteration.
        disp('   ');
        disp(['===== Current iteration number: ',num2str(obj.iiter),' =====' ]);
        disp(['Objective value: ', num2str(obj.objval)]);
        if obj.plot_on
            drawnow;
            idx = ~isnan(obj.objvals_);
            plot(find(idx==1), obj.objvals_(idx));
            xlabel('Iteration Number');
            ylabel('Function value');
        end
    end  
end


methods (Access = private)
    function display_iters(obj)
    % Display result while solving the problem in iterations
        if ~mod(obj.iiter, obj.niter_display) == 0; return; end      
        obj.display_result();
    end

    function STOP = stop(obj)
    % Stop algorithm base on iteration number and criterias.
        STOP = false;
        if obj.iiter == obj.maxiter
            obj.status = 'Reaches maximum iteration number';
            STOP = true;
        elseif obj.stop_criterias()
            obj.status = 'Meets all stopping criteria';
            STOP = true;
        end
    end
end
    
end % classdef