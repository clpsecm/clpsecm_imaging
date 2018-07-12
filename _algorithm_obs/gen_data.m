%% ===== Parameter Definition ===== %%
m = [64,64]; %-sparse image dimension
k = [8,8]; %-kernel size (Gaussian)
n = m+k-1; %-observed matrix size
theta = 0.002; %-sparsity level
nangles = 5;
XWorldLimits = [-1,1];
YWorldLimits = [-1,1];


%-Construct Gaussian kernel D
[X1,X2] = meshgrid(linspace(-1,1,k(1)),linspace(-1,1,k(2)));
Sigma = 0.3*eye(2);
Mu = [0,0];
D = mvnpdf([X1(:),X2(:)],Mu,Sigma);
D = reshape(D,size(X1));

%-Construct sparse matrix X
Xstar = zeros(m);
for I = 1:numel(Xstar)
    Xstar(I) = rand()<theta;
end

%-Construct image Y
Y = conv2(D,Xstar);
RefY = imref2d(size(Y),XWorldLimits,YWorldLimits);


%-Construct Observed image Y with random angles
angles = pi*rand(nangles,1)-pi/2;
% RY = line_integral(Y,angles);
RY = fourier_line_integral(Y,angles);

q = 1:nangles;
% Yhat = line_integral_adjoint(RY(:,q),angles(q),n);
LtRY = fourier_line_integral_adjoint(RY(:,q),angles(q),n);

%-Plot data
figure(1);
subplot(221); spy(Xstar); title('X')
subplot(222); imagesc(Y); title('Y')
subplot(223); imagesc(RY); title('RY')
subplot(224); imagesc(LtRY); title('L^*[RY]');

