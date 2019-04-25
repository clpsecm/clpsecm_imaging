% EXAMPLE3 Demostrates how to generate a simulated SECM data, and perform
% image reconstruction with reweighted ipalm using the CLP data from file
% clpconfig.  
close all;
initpkg

%% ===== Parameter Setting ====== %%
ticks  = [-3:0.05:3]; % Distance tags of each line measurment
ndiscs = 4;           % Number of discs
d_rad  = 0.1;         % Disc radius (mm)


%% ===== Generate discs and ground truth image ====== %%
% Generate disc D
func = @(x,y) x.^2 + y.^2 < d_rad^2;
D = DictProfile(ticks, func);

% Generate random map X0
X0 = SparseMap(ticks, 'random-location', d_rad, ndiscs);

% Generate simulated image Y
Y = X0 * D; 

% Generate lines R from image Y using clpconfig parameter
R = Y.line_project();

%% ===== Reconstrct the image with ReweightIPalm algorithm package ====== %%
% Generate problem with formulation 'Calibrated Lasso'
lda = 0.2*max(pos(D*back_project(R)));
prb = CalibLasso(R,D,R.params,lda);

% Generate algorithm 'ReweightIPalm' and solve the simulated problem.
alg = ReweightIPalmSecmSimul(prb,X0,D);
figure; 
alg.solve();

 
