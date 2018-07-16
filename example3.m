% EXAMPLE3 Demostrates how to generate a simulated SECM data, and perform
% image reconstruction when parameters are unkonwn.
close all;
initpkg

%% ===== Parameter Setting ====== %%
ticks = [-1:0.01:1];  % Distance tags of each line measurment
ndiscs = 6;           % Number of discs
d_rad = 0.05;          % Disc radius (mm)


%% ===== Generate discs and ground truth image ====== %%
% Generate disc D
func = @(x,y) x.^2 + y.^2 < d_rad^2;
D = DictProfile(ticks, func);

% Generate random map X0
X0 = SparseMap(ticks, 'random', d_rad, ndiscs);

% Generate simulated image Y
Y = X0 * D;

% Generate lines R from image Y using clpconfig parameter
R = Y.line_project();

%% ===== Reconstrct the image with IPalm algorithm package ====== %%
% Assign a WRONG parameter
p_wrong = ProbeParams('config',ticks);
p_wrong = p_wrong + {'angles',1};

% Solve the Calibration problem, start from wrong parameters
prb = CalibLasso(R,D,p_wrong);
alg = IPalmSecmSimul(prb,X0,D);
figure; 
alg.solve();

 
