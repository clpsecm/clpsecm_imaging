function TdtY = fftconv2_adjoint(Y,D)

[ky,kx] = size(D);
[ny,nx] = size(Y);
D0 = [D zeros(ky,nx-kx) ; ...
      zeros(ny-ky,nx) ];

TdtY = ifft2( fft2(Y) .* conj(fft2(D0)) );
TdtY = TdtY(1:ny-ky+1,1:nx-kx+1);