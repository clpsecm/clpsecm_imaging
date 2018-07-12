function RY = fourier_line_integral(Y,angles,varargin)
% [Description]: Discrete line integral by subpixel projection
%
% [Input] Y(ny,nx):         input image
%         angles(nangles):  in [-pi,pi); angle of lines
%         nbins:            number of bins after projection, 
%         offset(nangles):  projection of rotational origin to lines
%         Cy(2,1):          [y,x]; center of Y relative to rotational origin
%
% [Output]: RY(nbin,nangles) projected output of image
%
% [Algorithm]: 
%      L_t[Y] = L_0[R_{-t}[Y]] = iF{ S_0[ R_{-t}[F{Y}] ]}
%   thus we apply inverse rotation to fourier image FY
%
%   Variable offset is the relative location of projected rotational center
%   onto a line w.r.t. center of line
%
%   To rotate, we apply rotation by shearing:
%      [ cos(t) -sin(t) ] = [    1    0] * [ 1 -sin(t)] * [   1     0]
%      [ sin(t)  cos(t) ]   [tan(t/2) 1]   [ 0    1   ]   [tan(t/2) 1]
%             Rt          =  S(tan(t/2)) *  S'(-sin(t)) *  S(tan(t/2)) 
%   where S represents vertical shearing (y-axis) transform
%
% [Important Note]:
%  *** Relative axis of image ***
%   When rotation center is (0,0), boundary of image is: 
%        [-N/2-1/2, N/2-1/2] 
%   When line has offset 0, boundary of line is:
%        [-fl(nbins/2)-1/2,cl(nbin/2)-1/2]




%% ===== Parameter Defintion ====== %%
%-Zero padding Y
%-zeropad Y N-sqrare, where N = 2*max(size(Y)) 
[ny,nx] = size(Y);
N = 2*max(ny,nx);
Y0 = zeros(N,N);
Y0(N/2-floor((ny-1)/2):N/2+ceil((ny-1)/2),...
   N/2-floor((nx-1)/2):N/2+ceil((nx-1)/2) ) = Y;
Y = Y0;


%-nangles
nangles = length(angles);


%-nbins, centers, centerY
if nargin == 2
    nbins = N;
    offset = zeros(nangles,1);
    Cy = [0;0];
elseif nargin == 3
    nbins = varargin{1};
    offset = zeros(nangles,1);
    Cy = [0;0];
elseif nargin == 5
    nbins = varargin{1};
    offset = varargin{2};
    Cy = varargin{3};
end


%-Output RY
RY = zeros(nbins,nangles);


%---Fourier Rotation Opeartor---%
%-Rotation Center at pixel [N/2+1,N/2+1]
M = [0:N/2-1 0 -N/2+1:-1]';
shearv = @(Y,m)real( ifft( fft(Y) .* exp(-m*2i*pi*(M*M')/N) ) );
shearv = @(Y,m)fftshift(shearv(fftshift(Y),m));
shearh = @(Y,m)shearv(Y',m)';
rotation = @(Y,t)shearv( shearh( shearv(Y,tan(t/2)) ,-sin(t)) ,tan(t/2));
rotation = @(Y,t)rotation(rotation(Y,t/2),t/2); %-Smaller angle less error
rotation = @(Y,t)rotation(rotation(Y,t/2),t/2); %-Smaller angle less error




%% ===== Line Integral via Rotation by Shearing ===== %%
% close all;
for I = 1:nangles
    %-Line integral
    t = angles(I);
    L = sum(rotation(Y,t));
    LCy = cos(t)*Cy(2) + sin(t)*Cy(1);
    
    %-Calculate absolute location of line from centerY/offset
    %-     |-------(|c|)-------|      <-L; length N     
    %    |---(Off)--------|           <-RY; length nbins  
    %  <-{-(---0----------}----)--->  <-X axis
    leftBdryL   = round(LCy - N/2)-0.5; 
    rightBdryL  = leftBdryL + N;
    leftBdryRY  = round(-offset(I)-floor(nbins/2))-0.5;
    rightBdryRY = leftBdryRY + nbins;
    
     
    if leftBdryL <= leftBdryRY && rightBdryRY <= rightBdryL
        %-RY contained in L
        S = leftBdryRY - leftBdryL; 
        RY(:,I) = L(S+1:S+nbins);
    elseif leftBdryL <= leftBdryRY && rightBdryL < rightBdryRY
        %-RY outside of L on right side
        S = leftBdryRY-leftBdryL; 
        RY(1:N-S,I) = L(S+1:N);
    elseif leftBdryRY < leftBdryL && rightBdryRY <= rightBdryL 
        %- RY outside of L on left side
        S = leftBdryL - leftBdryRY;
        RY(S+1:nbins,I) = L(1:nbins-S);
    elseif leftBdryRY < leftBdryL && rightBdryL < rightBdryRY
        %- RY outside of L on both side
        S = leftBdryL - leftBdryRY;
        RY(S+1:S+N,I) = L;
    end

end
RY = RY/(nbins);

    
