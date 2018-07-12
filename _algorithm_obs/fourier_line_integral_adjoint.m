function Y = fourier_line_integral_adjoint(RY,angles,n,varargin)
% [Description]: Discrete line integral by subpixel projection
%
% [Input] RY(nbins,nangle): lines
%         angles(nangles):  vectors in [-pi,pi), angle of lines
%         n(2,1):           [ny,nx], dimension of output Y 
%         offset(nangles):  projection of rotational origin to lines
%         Cy(2,1):          [y,x]; center of Y relative to rotational origin
%
% [Output]: Y(ny,nx): output image 
%
% [Algorithm]: 
%      L_t[Y] = L_0[R_{-t}[Y]] = iF{ S_0[ R_{-t}[F{Y}] ]}
%   thus we apply inverse rotation to fourier image FY
%
%   To rotate, we apply rotation by shearing:
%      [ cos(t) -sin(t) ] = [    1    0] * [ 1 -sin(t)] * [   1     0]
%      [ sin(t)  cos(t) ]   [tan(t/2) 1]   [ 0    1   ]   [tan(t/2) 1]
%             Rt         =  S(tan(t/2)) *  S'(-sin(t)) *  S(tan(t/2)) 
%   where S represents vertical shearing (y-axis) transform


%% ===== Parameter Defintion ====== %%

%-Input Management
ny = n(1);
nx = n(2);

%-nbins,nangles, N
[nbins,nangles] = size(RY);
N = 2*max(n);

%-offset,centerY
if nargin == 3
    offset = zeros(nangles,1);
    Cy = [0;0];
elseif nargin == 5
    offset = varargin{1};
    Cy = varargin{2};
else
    error('Incorrect Number of Input Variables');
end

%-Output RY
Y = zeros(n);

%---Fourier Rotation Opeartor---%
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
    t = angles(I);
    %-Calculate absolute location of line from centerY/offset
    %-     |-------(|c|)-------|      <-L; length N     
    %    |---(Off)--------|           <-RY; length nbins  
    %  <-{-(---0----------}----)--->  <-X axis
    LCy = cos(t)*Cy(2) + sin(t)*Cy(1);
    leftBdryL   = round(LCy - N/2)-0.5; 
    rightBdryL  = leftBdryL + N;
    leftBdryRY  = round(-offset(I)-floor(nbins/2))-0.5;
    rightBdryRY = leftBdryRY + nbins;
    
    L = zeros(N,1);
    if     leftBdryL <= leftBdryRY && rightBdryRY <= rightBdryL
        %-RY contained in L
        S = leftBdryRY - leftBdryL; 
        L(S+1:S+nbins) = RY(:,I);
    elseif leftBdryL <= leftBdryRY && rightBdryL  < rightBdryRY
        %-RY outside of L on right side
        S = leftBdryRY-leftBdryL; 
        L(S+1:N) = RY(1:N-S,I);
    elseif leftBdryRY < leftBdryL  && rightBdryRY <= rightBdryL 
        %- RY outside of L on left side
        S = leftBdryL - leftBdryRY;
        L(1:nbins-S) = RY(S+1:nbins,I);
    elseif leftBdryRY < leftBdryL  && rightBdryL  < rightBdryRY
        %- RY outside of L on both side
        S = leftBdryL - leftBdryRY;
        L = RY(S+1:S+N,I);
    end
    
    %-Adjoint of Line Integral
    YI = ones(N,1)*L';
    YI = rotation(YI,-angles(I));
    Y = Y + YI(N/2-floor((ny-1)/2):N/2+ceil((ny-1)/2),...
               N/2-floor((nx-1)/2):N/2+ceil((nx-1)/2));
end
Y = Y/(nbins);

    
