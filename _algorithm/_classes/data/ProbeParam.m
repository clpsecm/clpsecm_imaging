classdef ProbeParam < handle
% PROBEPARAM Constrct a single parameter object for line probe.
% 
% obj = PROBEPARAM(value,bound,func) Given a 'value' vector(nval,1), the  
% 'bound' matrix is be either identical to value vector (constant value),
% or a matrix [LB of value , UB of value] (variable value). Input 'func'  
% should be a function takes in value vector.
% 
% obj = PROBEPARAM(NaN) Construct inactive PROBEPARAM
%     = PROBEPARAM(NaN,NAN,@(v)v)
%
% obj = PROBEPARAM(value) Construct constant value PROBEPARAM 
%     = PROBEPARAM(value,value,@(v)v)
%
% obj = PROBEPARAM(obj1) Copy PROBEPARAM.
%
% Ex:  Create PROBEPARAM with varying column magnitude
%   obj = PROBEPARAM( [1;2;3], [[0;1;2],[2;3;4]], @(v) ones(10,1)*v' ) 
% 
% PROBEPARAM is a handle object.
% 
% PROBEPARAM methods:
%   PLUS        -  Overload '+', Sum value of two input object.
%   MINUS       -  Overload '-', Subtract value between two input object.
%   MTIMES      -  Overload '*', Scalar multiple of value.
%   NORM        -  Return 2-norm of value.
%   INPROD      -  Return inner product between values.
%   SET_VALUE   -  Set value for parameter.
%   SET_BOUND   -  Set bound of value with input matrix or function handle.
%   SET_FUNC    -  Set function for parameter.
%   BOUND_VALUE -  Bound the value w.r.t. bound. 
%
% PROBEPARAM public fields:
%   VALUE       -  The value of parameter
%   BOUND       -  The upper/lower limit of all values
%   FUNC        -  A function handle produces function from value
%   ISVARIABLE  -  True if value is variable. False if value is constant.
%   ISACTIVE    -  True if the parametric operation on CLP is activated.
%
% Notice that obj.isvaraible is true iff obj.value equals obj.bound. 
%             obj.isactive   is true iff obj.value is not NaN.
%
% See also PROBEPARAMS

properties (SetAccess = private)    
    isvariable % Logical; True if value ~= bound
    isactive   % Logical; True if value ~= NaN
    value      % vector(n,1) The value of parameter
    bound      % vector(n,2) The upper/lower limit of all values
    func       % Function handle, produces parametric function from value.
end

methods
    function obj = ProbeParam(varargin)
    % obj = PROBEPARAM(value,bound,func) Construct parameter object
    % obj = PROBEPARAM(obj1)  Copy contructor 
    % obj = PROBEPARAM(value) Constant constructor
    % obj = PROBEPARAM(NaN)   Non-active constructor
        if nargin == 1
            v = varargin{1};
            if isa(v,'ProbeParam')
                obj.value = v.value;
                obj.bound = v.bound;
                obj.func = v.func;
                obj.isvariable = v.isvariable;
                obj.isactive = v.isactive;
            elseif isnan(v) 
                obj = ProbeParam(NaN,NaN,@(v)v);
            else % input contant v
                if isrow(v); v=v'; end
                obj = ProbeParam(v,v,@(v)v);
            end
        else % input (value,bound,func)
            obj.set_value(varargin{1});
            obj.set_bound(varargin{2});
            obj.set_func(varargin{3});
        end
    end

    function obj = plus(obj1,obj2)
    % obj = ProbeParam + ProbeParam  Return obj with summed values.
    % obj = ProbeParam + value       Return obj with summed values.
    % obj = value + ProbeParam       Return obj with summed values.
        if isa(obj1,'ProbeParam')
            obj = ProbeParam(obj1);
            if isa(obj2,'ProbeParam')
                check_issame(obj1,obj2);
                obj.set_value(obj.value + obj2.value);
            else % obj2 is value
                if isrow(obj2); obj2 = obj2'; end
                obj.set_value(obj.value + obj2);
            end
        else % obj1 is value
            if isrow(obj1); obj1 = obj1'; end
            obj = ProbeParam(obj2);
            obj.set_value(obj1 + obj.value);
        end
    end

    function obj = minus(obj1,obj2)
    % obj = ProbeParam - ProbeParam  Return obj with subtracted values.
    % obj = ProbeParam - value       Return obj with subtracted values.
    % obj = value - ProbeParam       Return obj with subtracted values.
        if isa(obj1,'ProbeParam')
            obj = ProbeParam(obj1);
            if isa(obj2,'ProbeParam')
                check_issame(obj1,obj2);
                obj.set_value(obj.value - obj2.value);
            else % obj2 is value
                if isrow(obj2); obj2 = obj2'; end
                obj.set_value(obj.value - obj2);
            end
        else % obj1 is value
            if isrow(obj1); obj1 = obj1'; end
            obj = ProbeParam(obj2);
            obj.set_value(obj1 - obj.value);
        end
    end

    function obj = mtimes(obj1,obj2)
    % obj = ProbeParam * scalar Returns obj with scalar multipled value.
    % obj = scalar * ProbeParam Returns obj with scalar multipled value. 
        if isa(obj1,'ProbeParam')
            obj = ProbeParam(obj1);
            obj.set_value(obj.value * obj2);
        else % obj2 is ProbeParam
            obj = ProbeParam(obj2);
            obj.set_value(obj.value * obj1);
        end
    end
    
    function n = norm(obj); n = norm(obj.value,2); end
    % n = obj.NORM() Return 2-norm of value
    
    function n = inprod(obj1,obj2)
    % n = obj.INPROD(obj1,obj2) Return inner product between values.
        check_issame(obj1,obj2);
        n = sum(obj1.value.*obj2.value);
    end
    
    function set_value(obj,value)
    % obj.SET_VALUE(value) Set value. isactive = true if value ~= NaN,
    % and isvariable = true if value ~= bound
        if isrow(value); value = value'; end
        if isnan(value)
            obj.isactive = false;
        else % value is not NaN
            obj.isactive = true;
        end
        obj.value = value;
        obj.check_isvariable;
    end
    
    function set_bound(obj,bound)
    % obj.SET_BOUND(bound) Set bound. Input bound can be a function handle
    % or matrix. isvariable = true if value ~= bound
        if isa(bound,'function_handle'); bound = bound(obj.value); end  
        if size(bound,1) ~= size(obj.value,1)
            error('param.bound should be of size (nval,1) or (nval,2)');
        end
        obj.bound = bound;
        obj.check_isvariable;
    end
    
    function set_func(obj,func); obj.func = func; end
    % obj.SET_FUNC(func) Set func.
    
    function bound_value(obj)
    % obj.BOUND_VALUE() Bound the value w.r.t. bound.
        if size(obj.bound,2) == 1
            obj.value = obj.bound;
        else % obj.bound is [LowerBound, UpperBound]
            v = obj.value;
            bl = obj.bound(:,1);
            bu = obj.bound(:,2);
            obj.set_value((bl>=v).*bl + (bu<=v).*bu + (bl<v & v<bu).*v);
        end
    end
end

methods (Access = private)
    function check_isvariable(obj)
        % Check if obj.value isvariable
        if isequal(obj.bound,obj.value) || isnan(obj.value(1))
            obj.isvariable = false;
        else
            obj.isvariable = true;
        end
    end
    
    function check_issame(obj1,obj2)
        if ~isequaln(obj1.bound,obj2.bound) || ...
           ~isequaln(obj1.func(obj1.value),obj2.func(obj1.value))
            error('Two parameter objects do not have same bound/func')
        end
    end
end

end % classdef

