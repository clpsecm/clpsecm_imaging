% EXAMPLE1 Reads data from file and introduce basic objectes in this code
% packages, including 

close all;
initpkg

%% ===== Read line data ===== $$
% Assign dataset
d = DataSpec('053018',0,4,1);    % (date, sample_number, nlines, version)
d.set_scan_distance(10,70,0.01); % (start_location, length, resolution)(mm)

d = DataSpec('071118',0,4,1);   % (date, sample_number, nlines, version)
d.set_scan_distance(30,40,0.01) 

%% ===== Probe parameter object ====== %%
% Set CLP parameters.
params        = ProbeParams(NaN);          % Initialize inactive parameter 
params.angles = ProbeParam([0,60,120,180]); % Set constant CLP angles.
params.shifts = ProbeParam([0,14.8,14.8,14.8]);

% Checkout all methods for class "ProbeParams"
whos params;      pause();
help ProbeParams; pause();


%% ===== Scan line object ====== %%
% Generate lines object, downsample 10 times.
% When input w/o parmas, read clpconfig for CLP parameter setting.
lines = d.get_clpsecm_data(params); 
lines.downsample(10,params);

% Checkout all methods for class "ScanLines"
whos lines;     pause();
help ScanLines; pause();


%% =====  SECM image object ===== %%
% Back project image from lines
bpimage = lines.back_project();

% Checkout all methods for class "SecmImage"
whos bpimage;   pause();
help SecmImage; pause();


%% ===== Plot data ===== %% 
figure();
subplot(121); lines.plot_lines();   title('Lines');
subplot(122); bpimage.draw_image(); title('Back project image');

