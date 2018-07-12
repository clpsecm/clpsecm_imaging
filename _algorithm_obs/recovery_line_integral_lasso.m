function Xhat = recovery_line_integral_lasso(RY,D,angles,n,lambda)

% [Description] 
%  Recovery of sparse chemical reactors in image domain
%  We observe the random projected convoluted signal as:
%         RY = L[D*X]   with   X sparse
%  Optimization problem: min lambda*|X|_1 +  1/2*|| RY - L[D*X] ||_F^2  -----(1)
%  Solve (1) with Primal forward-backward proximal approach:
%       X <- prox_{|.|_1/lambda}[(X-(LD)'(LDX-RY)]
%
% [Input]:
%   RY(nbins,nangles): Lines after integration
%   D:                 Kernel function              
%   angles(nangles,1): angles of line slopes
%   n(2,1):            [ny,nx]; size(Y)
%   lambda:            penalty variable
%
% [Output]:
%   Xhat(m,m):  Output sparsity pattern


if nargin < 4
    error('Insufficient input');
end


LINE_INTEGRAL_METHOD = 'Fourier';
niter = 500;
niterDisplayImage = 5;
niterDisplayDetail = 350;


%% ===== Parameter Definition ===== %%
%-Dimension Variables
m = n-size(D)+[1,1];
nbins = size(RY,1);


%-Linear Operators and objective function
if strcmp(LINE_INTEGRAL_METHOD,'Direct')
    LD =  @(X,D,angles) line_integral(conv2(X,D),angles,nbins);
    LDt = @(RY,D,angles,n) conv2_adjoint(line_integral_adjoint(RY,angles,n),D);
elseif strcmp(LINE_INTEGRAL_METHOD,'Fourier')
    LD =  @(X,D,angles) fourier_line_integral(conv2(X,D),angles,nbins);
    LDt = @(RY,D,angles,n) conv2_adjoint(fourier_line_integral_adjoint(RY,angles,n),D);
end
f = @(RY,X,LDX,lambda)  lambda*sum(sum(abs(X))) + ...
       1/2*norm(RY-LDX,'fro')^2;


%-Lipschtz of gradient of smooth function
Lip = 20*norm(LD(D,D,angles),'fro').^2;
t = 1/Lip;

%-Allocate lambda
if nargin < 5
    lambda = .4*max(max(abs(LDt(RY,D,angles,n))));
end

%-Adjoint of line integrals
LtRY = LDt(RY,D,angles,n);


%% ===== Fista ===== %%
Xhatprev = zeros(m);
Yhat = zeros(m);
kprev = 1;


for J = 1:niter
    Xhat = Yhat - t*LDt( LD(Yhat,D,angles)-RY ,D,angles,n ); %-Forward Gradient
    Xhat = sign(Xhat).*max(abs(Xhat)-lambda*t,0); %-Backward Proximal Step
    k = (1+sqrt(1+4*kprev^2))/2;
    Yhat = Xhat + (kprev-1)/k*(Xhat-Xhatprev);

    Xhatprev = Xhat;
    kprev = k;

    if mod(J,niterDisplayImage) == 0
        fig = figure(2);
        fig.WindowStyle = 'docked';
        drawnow;
        subplot(221); imagesc(LtRY); title('L^t[RY]'); axis equal; 
        subplot(222); imagesc(Xhat); title('X_{est}'); axis equal;
        subplot(223); imagesc(conv2(Xhat,D)); title('Y_{hat}'); axis equal;
        % subplot(224); plot(D(round(size(D,1)/2),:)), title('D Slice');
        subplot(224); spy(abs(Xhat)>0.2); title('supp[X_{hat}]');
    end
    disp(['===== Number of Iteration : ', num2str(J), ' ====='])
    if mod(J, niterDisplayDetail) == 0;
        LDX = LD(Xhat,D,angles);
        disp(['f = ', num2str( f(RY,Xhat,LDX,lambda) ) ]);
        disp(['|RY-L[D*Xhat]| = ', num2str(norm(RY-LDX,'fro')) ]);
    end
end


