function [c,idx] = closest(array,x)
% [c,idx] = CLOSEST(array,x) Find closest value c and its index idx in 
% array from x.
[~,idx] = min(abs(array-x));
c = array(idx);
end

