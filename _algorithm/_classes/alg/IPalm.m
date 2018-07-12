classdef IPalm < Solver
% IPALM Solves the following problem 
%        min_{x1,..xn}  sum_i^n fi(xi) + h(x1,...,xn)
%   where 
%      1. fi are regulation convex functions of variable xi.
%      2. h is a smooth coupling function.
%
%   At k-th iteration, IPALM solver [1] operates:
%   For i = 1 to n:
%      accelerate xi:  yi_k     <- xi_k + alpha*( xi_k - xi_(k-1) )
%      partial gradhi: dhi      <- D_{xi}[h](yi_k, x(n\i)_k)
%      prox of fi:     xi_(k+1) <- prox_{fi,ti_k}[ yi_k - ti * dhi ] 
%   
%   uisng the backtracking stepsize [2] at proximal gradient step.  
%
% obj = IPALM(problem) Construct a solver for input problem. The problem
% should contains properties {funcf, funch, proxf, gradh}.
%
% IPALM methods:
%   SET_ALPHA        - Set acceleration parameter alpha.
%   SET_NITER_OBJVAL - Set number of iteration to calculate objval.
%   SET_BACKTRACK    - Set backtracking parameters
% 
% IPALM public fields:
%   FUNCF - vars->scalar; regulation convex functions.
%   FUNCG - vars->scalar; smooth coupling function.
%   PROXF - cell(nvars,1); proxf{I}: (vi,ti)|-> vi; prox of fi.
%   GRADH - cell(nvars,1); gradh{I}: v  |-> (dh,t); gradient/stepsize of h.
%  
% For other methods/fields, see object Solver.
%
% See also SOLVER
%
% References:
% [1] P. Thomas, S. Sabach. "Inertial Proximal Alternating Linearized
%     Minimization (iPALM) for nonconvex and nonsmooth problems."
% [2] A. Beck, M. Teboulle, "A fast iterative shrinkage-thresholding 
%     algorithm for linear inverse problems." 

properties 
    funcf % vars->scalar; regulation function f
    funch % vars->scalar; smooth coupling function h 
    proxf % cell(2,1); proxf{I}: (vi,ti)|-> vi; prox of fi.
    gradh % cell(2,1); gradh{I}:    v   |-> (dh,t); gradient/stepsize of h.  
end

properties (Hidden)
    vars_  % cell(nvars,1); variable of optimization, last iterate
    grads  % cell(nvars,1); gradient of vars, current iterate
    alpha  % scalar; momentum constant
    tmax   % scalar; back-tracking maximum stepsize.
    tdiv   % scalar; back-tracking stepsize shrink ratio.
    niter_objval; % scalar; niter to set objective value
end

methods
    function obj = IPalm(problem)
    % Constructor of IPLAM solver
        obj = obj@Solver(problem.vars);
        obj.vars_  = problem.vars;
        obj.funcf  = problem.funcf;
        obj.funch  = problem.funch;
        obj.gradh  = problem.gradh;
        obj.proxf  = problem.proxf;
        obj.objval = obj.funcf(obj.vars)+obj.funch(obj.vars);  
        obj.grads  = cell(obj.nvars,1);       
        obj.alpha  = 0.9;
        obj.tmax   = 100;
        obj.tdiv   = 2;
        obj.niter_objval = 20;
    end

    function set_niter_objval(obj,niter); obj.niter_objval = niter; end
    % obj.SET_MITER_OBJVAL(niter) Set number of iteration to update objval.
    
    function set_backtrack(obj,tmax,tdiv) 
    % obj.SET_BACKTRACK(tmax,tdiv) Set backtracking parameters, the initial
    % step size start from (t_const)*(tmax), and shink by tdiv.
        obj.tmax = tmax; 
        obj.tdiv = tdiv; 
    end
    
    function set_alpha(obj,alpha)
    % obj.SET_ALPHA(alpha) Set momentum variable alpha.
        if alpha < 0 || alpha > 1; error('alpha should be in [0,1]'); end
        obj.alpha = alpha;
    end
end

methods (Access = protected)
    function iter(obj)
    % obj.ITER() Run a single iteration of IPALM
        v = cell(obj.nvars,1);
        for I = 1:obj.nvars
            % -------- Accelerated variables -------- %
            u_ = obj.vars;
            u_{I} = obj.vars{I} + obj.alpha*( obj.vars{I} - obj.vars_{I} );
            
            % --- Calculate gradient and stepsize --- %
            [dh,t_] = obj.gradh{I}(u_);
            obj.grads{I} = dh;            
            
            % ---- Backtracking proximal gradient ---- %
            u = obj.vars; 
            t = obj.tmax * t_;             
            while true
                % Calcluate Linearized proximal term
                u{I} = obj.proxf{I}( u_{I} - t*dh, t );
                % Determine the step size
                du = u{I} - u_{I};
                if t <= t_ 
                    v{I} = u{I}; break; % Constant step size             
                elseif obj.funch(u) <= obj.funch(u_) + inprod(du,dh) ...
                                        + 1/(2*t)*norm(du)^2 
                    v{I} = u{I}; break; % Backtrack step size
                end
                t = max( t/obj.tdiv, t_ ); 
            end
        end
        % Update all variables 
        obj.vars_ = deepcopy(obj.vars);
        obj.vars  = v;
        obj.set_objval();
    end

    function set_objval(obj)
    % Calculate objective value
        if mod(obj.iiter, obj.niter_objval) == 0
            obj.objval = obj.funcf(obj.vars) + obj.funch(obj.vars);
        end
    end

    function stopping = stop_criterias(obj)
    % Stopping criteria of ipalm, return stopping as logical.
        if mod(obj.iiter,10) == 0
            stopping = true;
            for I = 1:obj.nvars
                stopping = stopping && ...
                ( norm(obj.grads{I})/norm(obj.grads{I}) < 1e-3 ) && ...
                ( true );
            end
        else
            stopping = false; 
        end
    end
end


end % classdef



