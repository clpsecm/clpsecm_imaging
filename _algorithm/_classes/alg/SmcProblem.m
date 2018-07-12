classdef SmcProblem < handle
    % SMCPROBLEM is abstract of an optimization problem written as:
    %       min_{x1,...,xk}   sum_i^k fi(xi) + h(x1,x2,...,xk)
    % fi are convex regulation functions, h is smooth coupling function.
    %
    % SMCPROBLEM is a handle object.
    
properties (Abstract)
    vars  % cell(k,1); variables to optimize.
    funcf % vars->scalar; regulation function f
    funch % vars->scalar; smooth coupling function h 
    proxf % cell(2,1); proxf{I}: (vi,ti)|-> vi; prox of fi.
    gradh % cell(2,1); gradh{I}:    v   |-> (dh,t); gradient/stepsize of h.
end

end

