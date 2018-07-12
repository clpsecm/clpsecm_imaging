function line = get_line(varargin)
%%  Return ax+by+c = 0, line = [a,b,c]
%   Input: ('Point-Slope',p1,m)
%          ('Point,Point',p1,p2)

%----Input Allocation----%
METHOD = varargin{1}; % 'Point-Slope', 'Point-Point'

if strcmp(METHOD,'Point-Slope') 
    p1 = varargin{2};
    m = varargin{3};
    
elseif strcmp(METHOD,'Point-Point')
    p1 = varargin{2};
    p2 = varargin{3};
    m = (p2(2)-p1(2))/(p2(1)-p1(1)); %-dy/dx
else
    error('Incorrect METHOD');
end

%---Find line with point-slope---%
if    (m == inf) || (m == -inf), a=1; b=0;
elseif m == 0,                   a=0; b=1;
else
    a = sign(m)*m; b = -sign(m);
    n = norm([a,b]);
    a=a/n; b=b/n;
end
c = -(a*p1(1) + b*p1(2));

%---Return line---%
line = [a,b,c];


