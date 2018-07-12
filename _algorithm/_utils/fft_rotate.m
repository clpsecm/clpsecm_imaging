function rotated_img = fft_rotate(img,angle)
% rotated_img = FFT_ROTATE(img,angle) Rotate imput squre image with angle 
% (deg) in counterclockwise direction. The rotation center is image center.
% 
% To rotate, we apply rotation by shearing [1]:
%   [ cos(t) -sin(t) ] = [    1    0] * [ 1 -sin(t)] * [   1     0]
%   [ sin(t)  cos(t) ]   [tan(t/2) 1]   [ 0    1   ]   [tan(t/2) 1]
%           Rt         =  S(tan(t/2)) *  S'(-sin(t)) *  S(tan(t/2))
%
% Shear transform can be done via fft. The y-shear (column) transform 
% I(S_y(m)) has:
%        I(S_y(m)) =  iF_y{ F_y{I}(x,u) .* exp(-i2pi*m*x*u) }
%
% [1] K. Larkin et. al. Fast Fourier method for the accurate rotation of 
%     sampled images. 

[M,N] = size(img);
if M~=N; error('Input matrix should be square.'); end

%----Restrict angle between [-45º,45º] (1st orth)----%
angle = mod(angle+45,360)-45;
if     (angle>=-45) && (angle<= 45);                    orth = 1;
elseif (angle > 45) && (angle<=135); angle = angle-90;  orth = 2; 
elseif (angle >135) && (angle< 225); angle = angle-180; orth = 3;
elseif (angle>=225) && (angle< 315); angle = angle-270; orth = 4;                        
end

%----Roatete by shearing----%
k = ifftshift( [-floor(N/2) : (floor(N/2)+mod(N,2)-1) ]' );
shearv = @(Y,m) real(fftshift(ifft(fft(ifftshift(Y)).*exp(-m*(2i*pi)*(k*k')/N))));
shearh = @(Y,m) shearv(Y',m)';
rotation = @(Y,t) shearv( shearh( shearv(Y,tan(t/2)) ,-sin(t)) ,tan(t/2));
 
%----Rotate by angle + 90*orth----%
t = angle*(pi/180);
switch orth
    case 1; rotated_img = rotation(img,t);
    case 2; rotated_img = rot90(rotation(img,t),3); % rotate 90º
    case 3; rotated_img = rot90(rotation(img,t),2); % rotate 180º
    case 4; rotated_img = rotation(rot90(img,1),t); % rotate -90º
end
  