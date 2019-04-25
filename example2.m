% EXAMPLE2 Demostrates how to generate a simulated SECM data, and perform
% image reconstruction when all parameters are known.
close all;
initpkg

%% ===== Parameter Setting ====== %%
ticks = [-1:0.01:1];   % Distance tags of each line measurment (mm)
ndiscs = 30;           % Number of discs
disc_radius = 0.05;    % Disc radius (mm)
angles = [0:20:160]';  % Scan angles 


%% ===== Generate discs and ground truth image ====== %%
% Parameter of data, enable angles parameter only.
p        = ProbeParams(0);     % Initialize 0 parameters.
p.angles = ProbeParam(angles); % Set constant parameter angles.

% Generate disc D
func = @(x,y) x.^2 + y.^2 < disc_radius^2;
D = DictProfile(ticks, func);

% Generate random map X0
X0 = SparseMap(ticks, 'random-location', disc_radius, ndiscs);

% Generate simulated image Y
Y = X0 * D;

% Generate lines R from image Y
R = Y.line_project(p);

% Generate back projection image LtR from lines R
LtR = R.back_project();

% Plot results
figure;
subplot(221); Y.draw_image();          title('True Image');
subplot(222); D.draw_image();          title('Kernel');
subplot(223); R.plot_lines([0,40,80]); title('Lines');
subplot(224); LtR.draw_image();        title('Back projection image');


%% ===== Reconstrct the image with IPalm algorithm package ====== %%
% Generate problem with formulation 'Calibrated Lasso'
lda = 0.5*max(pos(D*D));
prb = CalibLasso(R,D,p,lda);

% Generate algorithm 'IPalm' and solve the simulated problem.
alg = IPalmSecmSimul(prb,X0,D);
figure; 
alg.solve(); 
 
