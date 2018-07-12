function [TdX] = fftconv2(X,D)

[ky,kx] = size(D);
[my,mx] = size(X);
X0 = [X zeros(my,kx-1) ; ...
       zeros(ky-1,mx+kx-1) ];
D0 = [D zeros(ky,mx-1) ; ...
       zeros(my-1,kx+mx-1) ];

TdX = ifft2(fft2(X0).*fft2(D0));





