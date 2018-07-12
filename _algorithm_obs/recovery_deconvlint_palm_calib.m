function [Xhat,D,params] = recovery_deconvlint_palm_calib...
         (RY,params_init,dim,lambda) 
% [Description] 
%  Recovery of sparse chemical reactors in image domain
%  We observe the random projected convoluted signal as:
%     RY = L[D*X]   with   X sparse
%  Optimization problem: 
%     min_{X,p} lambda*|X|_1 + 1/2*|| RY - L_p[D*X] ||_F^2  -----(1)
%  Solve (1) with alternating proximal gradient descent
%
%  a. X   <- prox_{|.|_1/(t1*lambda)}[ X - t1*(LD)'(LDX - RYtau) ]
%  b. tau <- tau - t2*Jtau'(RYtau - LDX)
%     theta <- theta - t3*Jtheta'(LDX-RYtau)
%     rad <- rad - t4*Jrad'(LDX-RYtau)
%     sigma <- sigma - t5*Jsigma'*(LDX - RYtau)
%     alpha <- alpha - t6*Jalpha'*(RYtau-LDX)
%
% [Input]:
%   RY(nbins,nangles): Lines after integration
%   LD: function handle 
%   params_init:  Initial angles of line slopes
%   dim:          Dimension parameter
%   lambda:       Penalty variable for lasso
%
% [Output]
%   Xhat(m):  Output sparse position
%   D      :  Kernel profile
%   params.amplitude:  amplitude for each lines in RY
%   params.rotation:   angle of slopes for each lines in RY
%   params.translation: transition for each lines in RY
%   params.diameter:   diameter for the kernel profile D
%   params.sigma:      smoothing parameter for kernel profile D

if nargin < 3
    error('Insufficient input.');
end

LINE_INTEGRAL_METHOD = 'Fourier';

%% ===== Input Parameters ===== %%
%-Dimension Variables
n = dim.n;
d = dim.d;
pixelLength = dim.PixelLength;
offset = dim.Offset;
Cy = dim.Cy;


nbins = size(RY,1);
nlines = length(offset);


%-Integration Method & Integral interval
dtau = 3; %-[NOTE]-Should be integer
dtheta = 0.01;
ddiam = 10;
dsigma = 10;
dalpha = 0.1;


%-Gradient Parameters
CALIBRATE_TRANSLATION = 1;
CALIBRATE_ROTATION = 1;
CALIBRATE_SCALING = 1;
CALIBRATE_RADIUS = 0;
CALIBRATE_VARIANCE = 0;


%-Iteration Parameters
niter = 1000;
niterDisplayImage = 10;
niterDisplayDetail = 1;

%-Variable Domain
tauRange = [-500,500]; %-tranlation (um)
thetaRange = [-15,15]; %-rotation (degree)
diamRange = [30,400];  %-diameter of disk (um)
sigmaRange = [10,200]; %-variance of gaussian 


tauRange = tauRange/pixelLength;
thetaRange = thetaRange*pi/180;

thetaInit = params_init.Rotation/180*pi;


%% ===== Parameter Definition ===== %%

%-Linear Operators and objective function
if strcmp(LINE_INTEGRAL_METHOD,'Direct')
    LD =  @(X,D,angles) line_integral(cropconv2(X,D),angles,nbins);
    LDt = @(RY,D,angles,n) cropconv2_adjoint(line_integral_adjoint(RY,angles,n),D);
elseif strcmp(LINE_INTEGRAL_METHOD,'Fourier')
    LD =  @(X,D,angles) fourier_line_integral(cropconv2(X,D),angles,nbins,offset,Cy);
    LDt = @(RY,D,angles,n) cropconv2_adjoint(fourier_line_integral_adjoint(RY,angles,n,offset,Cy),D);
end
f = @(RY,X,LDX,lambda)  lambda*sum(sum(abs(X))) + ...
       1/2*norm(RY - LDX,'fro')^2;


%-Variable Initilization
alpha = params_init.Amplitude;        %-Initial scaling
theta = params_init.Rotation/180*pi;  %-Initial angle of slopes
tau   = params_init.Translation;      %-Initial translation
diam  = params_init.Diameter;         %-Initial diameter
sigma = params_init.Sigma;



%-Profile D
gaudiam = ceil(4*sigma*2);
disk = fspecial('disk',diam/2);
gau  = fspecial('Gaussian',[gaudiam,gaudiam],sigma);
D = imresize(conv2(disk,gau), d ); %-Smooth & Resize the kernel
D = D/norm(D,'fro');


%-Lipschtz of gradient of smooth function
X1 = rand(n);
X1 = X1/norm(X1,'fro');
LD2 = norm(LDt(LD(X1,D,theta),D,theta,n));
t1 = .05/LD2; %-Stepsize for X

RY2 = norm(RY,'fro').^2;
t2 = .2/RY2; %-Stepsize for tau,theta

if mod(n,2)==0 %[W.A] Take care of even dimension n
    LX2 = 1e-5*norm(LDt(LD(X1,D,theta),X1(1:end-1,1:end-1),theta,d)); %-norm(Jdiam) ~ 1e-4
else
    LX2 = 1e-5*norm(LDt(LD(X1,D,theta),X1,theta,d)); %-norm(Jdiam) ~ 1e-4
end
t3 = 1/LX2;  %-Stepsize for diam, sigma


%-L1-penalty variable
if nargin ~= 4
    lambda = .5*max(max(abs(LDt(RY,D,thetaInit,n))));
end



%% ===== PALM ===== %%
Xhat = zeros(n);
LDXhat = zeros(size(RY));
RYtau = RY;

for I = 1:niter    
    %---Assign RYtau in new iteration---%
    if CALIBRATE_TRANSLATION
        for J = 1:nlines
            tt = tau(J);
            RYtau(:,J) = alpha(J)*imtranslate(RY(:,J),[0,tt]);
        end
    end
    
    if CALIBRATE_SCALING
        for J = 1:nlines
            a = alpha(J)*norm(RY(:,J),2)/norm(RYtau(:,J),2);
            RYtau(:,J) = a*RYtau(:,J);    
        end
    end
    
    
    %---Assign kernel D in new iteration---%
    gaudiam = ceil(4*sigma*2);
    disk = fspecial('disk',diam/2);
    gau  = fspecial('Gaussian',[gaudiam,gaudiam],sigma);
    D = imresize(conv2(disk,gau), d ); %-Smooth & Resize the kernel
    D = D/norm(D,'fro');

    
    if CALIBRATE_RADIUS        
        gaudiam = ceil(4*sigma*2);
        disk = fspecial('disk',(diam+ddiam)/2);
        gau  = fspecial('Gaussian',[gaudiam,gaudiam],sigma);
        Dddiam = imresize(conv2(disk,gau), d ); %-Smooth & Resize the kernel
        Dddiam = Dddiam/norm(Dddiam,'fro');
        
    end
    
    if CALIBRATE_VARIANCE
        gaudiam = ceil(4*(sigma+dsigma)*2);
        gau  = fspecial('Gaussian',[gaudiam,gaudiam],sigma+dsigma);
        Ddsigma = imresize(conv2(disk,gau), d ); %-Smooth & Resize the kernel
        Ddsigma = Ddsigma/norm(Ddsigma,'fro');       
    end
    
  
    %---Proximal Gradient for Xhat---%
    Xhat = Xhat - t1*LDt( LDXhat-RYtau,D,theta,n ); 
    Xhat = sign(Xhat).*max(abs(Xhat)-lambda*t1,0);
    LDXhat = LD(Xhat,D,theta);
    

    %---Jacobian---%
    if CALIBRATE_TRANSLATION 
        Jtau = [zeros(dtau,nlines) ; RYtau(1:end-dtau,:)] - RYtau;
    end
    if CALIBRATE_SCALING
        Jalpha = (RYtau*(1+dalpha)-RYtau)/dalpha;
    end
    if CALIBRATE_ROTATION
        Jtheta = (LD(Xhat,D,theta+dtheta)-LDXhat)/dtheta;
    end
    if CALIBRATE_RADIUS
        Jdiam = (Dddiam-D)/ddiam;
    end
    if CALIBRATE_VARIANCE
        Jsigma = (Ddsigma-D)/dsigma;
    end
    
    %---Gradient for tau, theta, diam, dsigma---%
    if CALIBRATE_TRANSLATION
        for J = 1:nlines
            tau(J) = tau(J) - t2*Jtau(:,J)'*( RYtau(:,J) - LDXhat(:,J) );
        end
    end
    if CALIBRATE_SCALING
        for J = 1:nlines
            alpha(J) = alpha(J) - t2*Jalpha(:,J)'*(RYtau(:,J) - LDXhat(:,J));
        end
    end
    if CALIBRATE_ROTATION
        for J = 1:nlines 
            theta(J) = theta(J) - t2/2*Jtheta(:,J)'*( LDXhat(:,J) - RYtau(:,J) );
        end 
    end 
    if CALIBRATE_RADIUS || CALIBRATE_VARIANCE
        LXtM  = LDt(LDXhat-RYtau,Xhat,theta,d);        
        if CALIBRATE_RADIUS
            diam  = diam  - t3*sum(sum(Jdiam.*LXtM));
        end
        if CALIBRATE_VARIANCE
            sigma = sigma - t3*sum(sum(Jsigma.*LXtM));
        end       
    end

    
    %---Projection to constraint set---%
    tau = max(tau,tauRange(1));
    tau = min(tau,tauRange(2));
    
    theta = max(theta,thetaInit+thetaRange(1));
    theta = min(theta,thetaInit+thetaRange(2));
    
    diam = max(diam,diamRange(1));
    diam = min(diam,diamRange(2));
    
    sigma = max(sigma,sigmaRange(1));
    sigma = min(sigma,sigmaRange(2));
    
    %-alpha in simplex * nangles
    alpha = max(alpha,0);
    alpha = min(alpha,nlines);
    alpha = alpha/sum(alpha)*nlines;
    
   
    if mod(I,niterDisplayImage) == 0
        if I == niterDisplayImage
            fig = figure;
            fig.WindowStyle = 'docked';
            figNumber = get(fig,'Number');
        end
        figure(figNumber);
        drawnow;        
        LtRYtau = fourier_line_integral_adjoint(RYtau,theta,n,offset,Cy);
        subplot(221); subplot(221); imagesc(LtRYtau); title('L^t[RY]'); axis equal;
        subplot(222); imagesc(Xhat); title('X_{est}'); axis equal;
        subplot(223); imagesc(cropconv2(Xhat,D)); title('Y_{hat}'); axis equal;
        subplot(224); plot(D(round(d(1)/2),:)); title('D shape');
    end
    disp(['===== Number of Iteration : ', num2str(I), ' ====='])
    if mod(I, niterDisplayDetail) == 0;
        disp(['tau(um) = ', num2str(tau*pixelLength)]);
        disp(['theta(deg) = ', num2str(theta*180/pi)]);
        disp(['alpha = ', num2str(alpha)]);
        disp(['f = ', num2str( f(RYtau,Xhat,LDXhat,lambda) ) ]);
        disp(['|RYtau-L[D*Xhat]| = ', num2str(norm(RYtau-LDXhat,'fro')) ]);
        disp(['diam(D) = ', num2str(diam),'  sigma = ', num2str(sigma)]);       
    end
    
end

params.Amplitude = alpha;
params.Rotation = theta*180/pi;
params.Translation = tau*pixelLength;
params.Diameter = diam;
params.Sigma = sigma;


