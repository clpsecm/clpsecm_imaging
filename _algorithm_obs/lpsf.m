function y = lpsf(p,x)
%-Convoluted Exponential
% c * exp(-a|t|^m) * (t+1)^n
%
% Input: p(1,3), [a,m,n]
%        p(1,4), [a,m,n,c]
%        p(1,5), [a,m,n,c,b]
%        x scalar, function domain


%-Input Parameter Allocation
xm = x(round(0.85*numel(x)/2)); %-Find 'mean' of x
if numel(p) == 5
    a=p(1); m=p(2); n=p(3); b=p(4); c=p(5);
elseif numel(p) == 4
    a=p(1); m=p(2); n=p(3); b=p(4); c=1;
elseif numel(p) == 3
    a=p(1); m=p(2); n=p(3); b=xm;   c=1;
else
    error('Invalid Parameter Dimension')
end


%-Dimension Parameter
dd = x(2) - x(1); %-Pixel Width
s = 0:dd:max(x); %-Integrating variable


%-Convolution with integral
y = zeros(size(x));
for I = 1:numel(x)
    y(I) = c*sum( exp(-a*abs(x(I)-b - s).^m).*((s+1).^n)  )*dd;
end


%-Scaling y
if numel(p) ~= 5
    y = y/(sum(y)*dd); % Integral of y = 1
end
    