function y = normalize(x,p)
% y = NORMALIZE(x,p) y = x / ||x||_p 
y = x/norm(x,p);
end

