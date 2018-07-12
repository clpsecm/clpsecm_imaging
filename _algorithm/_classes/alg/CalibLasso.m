classdef CalibLasso < SmcProblem
% CALIBLASSO Generates the sparse recovery problen for CLP-SECM using lasso
% objective and while calibrating the line measurments. 
%
%     min_{X,p} lda*sum(X+) + Ind{p in P} + 1/2*||L_p[D*X]-R||F^2 ---- (1)
%                      
% where Ind is the indicator function. X+ projects X to positive orthant.
% To solve (1), we implements the following functional properties:
%   
%     proxf{1}(X) = prox_{lda*sum(.+)}[X] = soft_lda[X+]
%     proxf{2}(p) = prox{Ind{.in P}}[p] = proj{p to P}
%     gradh{1}(X) = D*L_p'[ L_p[D*X] - R ]
%     gradh{2}(p) = Jp'*[ L_p[D*X-R] ]
%
% obj = CALIBLASSO(lines,dict,params) Initialize the problem by providing
% the observed lines, the dictionary profile, and the CLP parameter
% settings. It generates the objective (1) and its functional properties.
%
% CALIBLASSO is a handle object.
%
% CALIBLASSO public fields:
%   VARS  - Varables for optimization, includes the sparse map and params.
%   FUNCF - Regulation function in objective. 
%   FUNCH - Smooth coupling function in objective. 
%   PROXF - Proximal operator w.r.t. regulation function of f.
%   GRADH - Gradient of smooth coupling term h.
%   LDA   - The sparsity penalty variable.
%
% The implementation of fields should match algorithm requirement. Such as
% IPALM solver class.
%
% See Also SMCPROBLEM, IPALM
properties
    vars  % cell(2,1); initial values for variables to optimize.
    funcf % vars->scalar; regulation function f
    funch % vars->scalar; smooth coupling function h 
    proxf % cell(2,1); proxf{I}: (vi,ti)|-> vi; prox of fi.
    gradh % cell(2,1); gradh{I}:    v   |-> (dh,t); gradient/stepsize of h.
    lda   % scalar; sparsity regularizer
end

methods 
    function obj = CalibLasso(lines,dict,params)
        % Construct an instance of problem  
        D = dict;
        R = lines;
        L  = @(img,params) line_project(img,params);
        Lt = @(lines)      back_project(lines);
        
        % Calculate lda
        p0 = ProbeParams(params);
        p0.shifts = ProbeParam(NaN);
        p0.intensity = ProbeParam(NaN);
        lda = 0.5*max(pos(D*Lt(L(D,p0))));
        obj.lda = lda; 

        % Define variables (X,params)
        obj.vars = { SecmImage(R.ticks), params }; 
        
        % Defind objective function f,h
        obj.funcf = @(v) lda*sum(pos(v{1}));  
        obj.funch = @(v) 1/2*norm(L(D*v{1},v{2})-R)^2; 
        
        % Define derivative gradh and step size t
        D2 = sum(D)^2;
        L2 = R.nlines * R.nmeasures;
        obj.gradh{1} = @(v) deal( D*Lt(L(D*v{1},v{2})-R) , 1/(D2*L2) );
        obj.gradh{2} = @(v) obj.gradh_params(v,D,R); 
        
        % Define proximal operator proxf
        obj.proxf{1} = @(vi,ti) soft(pos(vi),ti*lda); 
        obj.proxf{2} = @(vi,ti) obj.proxf_params(vi); 
    end
end

methods (Access = private)
    function [dh,t] = gradh_params(obj,v,D,R)
        % Define linear maps
        L = @(img,params) line_project(img,params);
        Rdiff = L(D*v{1},v{2})-R;
        Rdiff = Rdiff.currents;
        properties = {'angles','shifts','intensity','psf'};
        
        % Initializae derivative dh and stepsize t
        dh = ProbeParams(v{2});
        for prop = properties
            param = prop{1};
            dh.set_value(param,0);
        end
        t = 1/eps;
        
        % Assign gradient for each variables
        dv = 0.1; % Used for calculate Jacobian numetically.
        for prop = properties
            param = prop{1};
            if any(strcmp(param,{'angles','shifts','intensity'}))
                % ----- 'Perline' variables ----- %
                if v{2}.(param).isvariable
                    % Jacobian of L[D*X,p] w.r.t. param
                    J = L(D*v{1},v{2}+{param,dv}) - L(D*v{1},v{2});
                    J = J.currents/dv;
                    % Derivative dh and stepsize t
                    dhv = sum(Rdiff .* J)';
                    t   = 1/(1/t + norm(J,'fro')^2);
                    dh.set_value(param,dhv);
                end
            elseif any(strcmp(param,'psf'))
                % --- 'Non-seperable' variables --- %
                if v{2}.(param).isvariable
                    nvar = numel(v{2}.(param).value);
                    dhv = zeros(nvar,1);
                    for I = 1:nvar
                        % Jacobian of L[D*X,p] w.r.t. param
                        dvi = zeros(nvar,1);
                        dvi(I) = dv;
                        J = L(D*v{1},v{2}+{param,dvi}) - L(D*v{1},v{2});
                        J = J.currents/dv;
                        % Derivative dh and stepsize t
                        dhv(I) = sum(sum(Rdiff .* J));
                        t      = 1/(1/t + norm(J,'fro')^2);
                    end
                    dh.set_value(param,dhv);
                end
            end
        end
    end % end function
    
    function vi = proxf_params(obj,vi)
        vi.bound_values();
    end
end

end % classdef