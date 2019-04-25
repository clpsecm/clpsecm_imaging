classdef ProbeParams < handle
% PROBEPARAMS Constrct the parameters object for line probe.
% 
% obj = PROBEPARAMS(NaN) Construct inactive PROBEPARAMS.
%
% obj = PROBEPARAMS(0) Construct zero PROBEPARAMS, only 'angles' is [0].
%
% obj = PROBEPARAMS(obj1) Copy PROBEPARAMS.
%
% obj = PROBEPARAMS('config',ticks) Constructor with clpconfig.
%
% PROBEPARAMS is a handle object.
%
% PROBEPARAMS methods:
%   PLUS           -  Overload '+', Sum values in all params
%   MINUS          -  Overload '-', Subtract values in all params.
%   MTIMES         -  Overload '*', Scalar multiply value in all params.
%   NORM           -  Return 2-norm of active parameters values.
%   INPROD         -  Return Inner product of active param values.
%   SET_VALUE      -  Set value of target ProbeParam
%   SET_BOUND      -  Set bound of target ProbeParam
%   SET_FUNC       -  Set func of target ProbeParam
%   GET_ISVARIABLE -  Return which parameters are variables
%   GET_ISACTIVE   -  Return which parameters are active
%   BOUND_VALUES   -  Bound all parameter values using its bound. 
%
% PROBEPARAMS public fields:
%   ANGLES        -  Scan angles of lines.
%   SHIFTS        -  Shifts of each lines.
%   INTENSITY     -  Intensity of each lines, the 'mask' of lines.
%   PSF           -  Point-spread-function of line probe. 
%
% See also PROBEPARAM, CLPCONFIG
properties
    angles     % obj ProbeParam, scan angles of lines 
    shifts     % obj ProbeParam, shifts of each lines.
    intensity  % obj ProbeParam, intensity of each lines.
    psf        % obj ProbeParam, point-spread-function of line probe.
end

methods
    function obj = ProbeParams(varargin)
    % CALIBPARAMS Construct an instance of this class
        plist = {'angles','shifts','intensity','psf'};
        if nargin == 1
            v = varargin{1};
            iscopy = isa(v,'ProbeParams');
            if     iscopy;   for p=plist; obj.(p{1})=ProbeParam(v.(p{1})); end
            elseif v == 0;   for p=plist; obj.(p{1})=ProbeParam(NaN);   end
            elseif isnan(v); for p=plist; obj.(p{1})=ProbeParam(NaN); end  
            else 
                error('Input should be either 0, NaN or ProbeParams.');
            end
            if v == 0; obj.angles = ProbeParam(0); end
            return;
        elseif strcmp(varargin{1},'config')
            ticks = varargin{2};
            % ----- Initial value ----- %
            for p = plist; obj.(p{1}) = ProbeParam(NaN); end
            % ----- Assign Parameters ----- %
            clpcfg = clpconfig(ticks);
            for prop = fields(clpcfg)'
                param = prop{1};
                value = clpcfg.(param).value;
                bound = clpcfg.(param).bound;
                func  = clpcfg.(param).func;
                obj.(param) = ProbeParam(value, bound, func);
            end
        end
    end
    
    function obj = plus(obj1,obj2)
    % obj = obj1 + obj2 Return sum of parameters.
    % obj = obj1 + {param1,value1,...} Return sum target values.
    % obj = {param1,value1,...} + obj2 Return sum target values. 
        if isa(obj1,'ProbeParams')
            obj = ProbeParams(obj1);
            if isa(obj2,'ProbeParams')
                for prop = properties(obj1)'
                    param = prop{1};
                    obj.(param) = obj1.(param) + obj2.(param);
                end
            else % obj2 = {param1,value1,param2,value2,...}
                for prop = properties(obj1)'
                    param = prop{1};
                    idx = find(strcmp(obj2,prop{1}));
                    if ~isempty(idx)
                        obj.(param) = obj.(param) + obj2{idx+1};
                    end
                end
            end
        else % obj1 = {param1,value1,param2,value2,...}
            obj = ProbeParams(obj2);
            for prop = properties(obj2)'
                param = prop{1};
                idx = find(strcmp(obj1,prop{1}));
                if ~isempty(idx)
                    obj.(param) = obj.(param) + obj1{idx+1};
                end
            end
        end  
    end

    function obj = minus(obj1,obj2)
    % obj = obj1 + obj2 Return subtraction of parameters.
    % obj = obj1 + {param1,value1,...} Return subtracted target values.
    % obj = {param1,value1,...} + obj2 Return subtracted target values. 
        if isa(obj1,'ProbeParams')
            obj = ProbeParams(obj1);
            if isa(obj2,'ProbeParams')
                for prop = properties(obj1)'
                    param = prop{1};
                    obj.(param) = obj1.(param) - obj2.(param);
                end
            else % obj2 = {param1,value1,param2,value2,...}
                for prop = properties(obj1)'
                    param = prop{1};
                    idx = find(strcmp(obj2,param));
                    if ~isempty(idx)                     
                        obj.(param) = obj.(param) - obj2{idx+1};
                    end
                end
            end
        else % obj1 = {param1,value1,param2,value2,...}
            obj = ProbeParams(obj2);
            for prop = properties(obj2)'
                param = prop{1};
                idx = find(strcmp(obj1,param));
                if ~isempty(idx)
                    obj.(param) = obj.(param) - obj1{idx+1};
                end
            end
        end  
    end   
    
    function obj = mtimes(obj1,obj2)
    % obj = obj1 * scalar  Return multiplied values of all parameters.
    % obj = scalar * obj2  Return multiplied values of all parameters.
        if ~isscalar(obj1) || ~isscalar(obj2) 
            error('ProbeParams allows only scale multiplication.');
        end
        if isa(obj1,'ProbeParams')
            obj = ProbeParams(obj1);
            for prop = properties(obj1)'
                param = prop{1};
                obj.(param) = obj1.(param) * obj2;
            end
        else % obj2 is ProbeParams
            obj = ProbeParams(obj2);
            for prop = properties(obj2)'
                param = prop{1};
                obj.(param) = obj2.(param) * obj1;
            end
        end
    end
    
    function n = norm(obj)
    % n = obj.NORM() 2-norm of active parameter values.
        n = 0;
        for prop = properties(obj)'
            param = prop{1};
            if obj.(param).isactive
                n = sqrt(n^2 + norm(obj.(param))^2) ;
            end
        end
    end
    
    function n = inprod(obj1,obj2)
    % n = obj.INPROD(obj1,obj2) inner product of active parameter values
        n = 0;
        for prop = properties(obj1)'
            param = prop{1};
            if obj1.(param).isactive && obj2.(param).isactive 
                n = n + inprod(obj1.(param),obj2.(param)); 
            end
        end
    end
    
    function set_value(obj,param,value); obj.(param).set_value(value); end 
    % SET_VALUE(obj,param,value) Set value for param
    
    function set_bound(obj,param,bound); obj.(param).set_bound(bound); end
    % SET_BOUND(obj,param,bound) Set bound for param
    
    function set_func(obj,param,func); obj.(param).set_func(func); end 
    % SET_FUNC(obj,param,func) Set func for param    
    
    function p = get_isvariable(obj)
    % p = obj.GET_ISVARIABLE() Get cell of variable parameter names 
        p = {};
        for prop = properties(obj)'
            if obj.(prop{1}).isvariable 
                p{end+1} = prop{1}; 
            end
        end
    end
    
    function p = get_isactive(obj)
    % p = obj.GET_ISACTIVE() Get cell of variable parameter names 
        p = {};
        for prop = properties(obj)'
            if obj.(prop{1}).isactive 
                p{end+1} = prop{1}; 
            end
        end
    end

    function bound_values(obj)
    % obj.BOUND_VALUES() Bound all parameter values using its bound. 
        for prop = properties(obj)'
            param = prop{1};
            if obj.(param).isvariable
                obj.(param).bound_value();
            end
        end
    end
end 
    
end % classdef


