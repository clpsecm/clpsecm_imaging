function TdtY = conv2_adjoint(Y,D)
% [Input]:  Y(ny,nx), image in covolved domain
%           D(ky,kx), kenel
% [Output]: TdtY(ny-ky+1,nx-kx+1)
% [Algorithm]: Cropped cross correlation 

[ky,kx] = size(D);
TdtY = xcorr2(Y,D);
TdtY = TdtY(ky:end-ky+1, kx:end-kx+1);





