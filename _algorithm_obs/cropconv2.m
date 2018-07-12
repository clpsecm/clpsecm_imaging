function Y = cropconv2(X,D)
% [Description]: Output Y has same dimension as X

    %-Dimension of D
    [dy,dx] = size(D);
    if (mod(dy,2)==0) || (mod(dx,2)==0)
        error('Dimension of D should be odd')
    end
    d2y = (dy-1)/2;
    d2x = (dx-1)/2;
    

    %-crop Y
    Y = conv2(D,X);
    Y = Y(d2y+1:end-d2y, d2x+1:end-d2x);
    
end