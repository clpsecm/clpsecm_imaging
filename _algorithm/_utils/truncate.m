function y = truncate(x,set)
%TRUNCATE_MID(x,set) truncate vector x, keep the set 
%   if set is vector, y = x(set)
%   if set is 'mid', y is center middle part of x
n = length(x);

if isnumeric(set)
    y = x(set);
elseif strcmp(set,'mid')
    if mod(n+1,2)~=0; error('no centered middle segment'); end
    y = x(floor((n-1)/4)+1:floor(3*(n-1)/4)+1);
end



