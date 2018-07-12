function TdtY = cropconv2_adjoint(Y,D)
% [Input]:  Y(ny,nx), image in covolved domain
%           D(ky,kx), kenel
% [Output]: TdtY(ny-ky+1,nx-kx+1)
% [Algorithm]: Cropped cross correlation 


%-Dimension of D
[dy,dx] = size(D);
if (mod(dy,2)==0) || (mod(dx,2)==0)
    error('Dimension of D should be odd')
end
d2y = (dy-1)/2;
d2x = (dx-1)/2;

%-cross correlation
CtY = zeros(size(Y) + size(D)- [1,1]);
CtY(d2y+1:end-d2y, d2x+1:end-d2x) = Y;
TdtY = xcorr2(CtY,D);
TdtY = TdtY(dy:end-dy+1, dx:end-dx+1);





