function vector = shift(vector,s)
% vector = SHIFT(vector,s) Shift vector by interger s indicies, 
if abs(s) >= length(vector)
    vector = zeros(size(vector));
    return;
end
if s > 0
    vector(1+s:end) = vector(1:end-s);
    vector(1:s) = 0;
else % s <= 0
    s = abs(s);
    vector(1:end-s) = vector(s+1:end);
    vector(end-s+1:end) = 0;
end

