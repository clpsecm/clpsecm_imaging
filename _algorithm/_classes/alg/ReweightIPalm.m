classdef ReweightIPalm < Solver
% REWEIGHTIPALM Solve multiple iterates of IPalm, each with reweighted lda.
% For iterates k = 2,3,... the lda is determined by 
%                   lda^(k) = 1./(X^(k-1)+eps)
% where X^(k-1) is the solution of sparse map in last iterate.
% 
% obj = REWEIGHTIPALM(problem) Construct a solver for input problem. The 
% problem should contains properties {funcf, funch, proxf, gradh}.
%
% REWEIGHTIPALM methods:
%   SET_NITER_IPALM - Set number of iteration in each run of IPalm.
%
% For other methods/fields, see object Solver.
% 
% See also IPALM CALIBLASSO SOLVER

properties(Hidden)
    f_X    % function handle; null input function return sparse vector.
    prb    % Object CalibLasso; construct REWEIGHTIPALM with the problem.
    ipalm  % Object Ipalm; algorithm in each run. 
    niter_ipalm  % scalar; niter in each run of ipalm.
end

methods
    function obj = ReweightIPalm(problem)
    % REWEIGHTIPALM Construct Reweighted IPalm solver
        obj = obj@Solver(problem.vars);
        obj.f_X = @() obj.vars{1}.image;
        obj.prb = problem;
        
        obj.plot_on = false; 
        obj.set_display(1);
        obj.set_maxiter(6);
        obj.set_niter_ipalm(50);
    end
    
    function set_niter_ipalm(obj,niter_ipalm) 
    % obj.SET_NITER_IPALM(niter_ipalm) Set niter in each run of ipalm.
        obj.niter_ipalm = niter_ipalm; 
    end
end

methods (Access = protected)
    function iter(obj)
    % obj.ITER() Run one instance of IPALM
        X = obj.f_X();
        if norm(X) ~= 0
            obj.prb.set_lda(1e-3*obj.objval./(X+eps));
            obj.prb.set_init(obj.vars);
        end
        obj.ipalm = IPalm(obj.prb);
        obj.ipalm.set_display(5);
        obj.ipalm.set_maxiter(obj.niter_ipalm);
        figure();
        obj.ipalm.solve();
        close;
        obj.vars = obj.ipalm.vars;
    end
    
    function set_objval(obj); obj.objval = obj.ipalm.objval; end

    function stopping = stop_criterias(obj); stopping = false; end
end

end

