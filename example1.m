% EXAMPLE1 Reads data from file and introduce basic objectes in this code
% packages, including 
close all;
initpkg

%% ===== Read line data ===== $$
% Assign dataset
d = DataSpec('100818',3,7,2);    
d.set_scan_distance(2.4, 4.2, 0.01);
angles = [0,45,70,90,115,135,180];
% 
% d = DataSpec('021319',32,9,1);  
% d.set_scan_distance(2.4, 5.0, 0.01);
% angles = [0,20,40,60,80,100,120,140,160];

%% ===== Probe parameter object ====== %%
% Set CLP parameters.
p        = ProbeParams(NaN);   % Initialize inactive parameter 
p.angles = ProbeParam(angles); % Set constant CLP angles.

%% ===== Scan line object ====== %%
% [Note] When input w/o parmas, read clpconfig for CLP parameter setting.
lines = d.get_clpsecm_data(p);
lines.downsample(4);
lines.params.angles = ProbeParam(angles); 
lines.params.psf = ProbeParam(NaN);


%% ===== SECM image object ===== %%
% Back project image from lines
bpimage = lines.back_project();


%% ===== Plot data ===== %% 
figure();
subplot(121); lines.plot_lines();   title('Lines');
subplot(122); bpimage.draw_image(); title('Back project image');

