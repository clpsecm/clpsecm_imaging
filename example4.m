% EXAMPLE4 Reads data from file and introduce basic objectes in this code
% packages, including 
close all;
initpkg

%% ===== Read line data ===== $$
% Assign dataset
% d = DataSpec('091718',12,5,1);
% d.set_scan_distance(2.4,7,0.01);
% angles = [0,45,90,135,180];

d = DataSpec('100818',3,7,2);
d.set_scan_distance(2.4, 4.2, 0.01);
angles = [0,45,70,90,115,135,180];
% 
% d = DataSpec('110718',0,7,1);
% d.set_scan_distance(2.4,6,0.01);
% angles = [0,45,70,90,115,135,180];

% d = DataSpec('021319',32,9,1);
% d.set_scan_distance(2.4, 4.8, 0.01);
% angles = [0,20,40,60,80,100,120,140,160];

 
%% ===== Scan line object ====== %%
% Read clpconfig for CLP parameter setting.
lines = d.get_clpsecm_data();
lines.downsample(4);
lines.zeroing();
lines.params.angles = ProbeParam(angles);
 
% Choose angle subset
% agl = [0,45,70,115,135];
% [~,ind] = ismember(agl,angles);
% lines.params.angles = ProbeParam(agl);
% lines = ScanLines(lines.ticks,lines.currents(:,ind),lines.params);


%% ===== Setup Signal Assumptions and Probe Parameters ===== %%
% -Setup disc D
disc_radius = 0.075;
func = @(x,y) x.^2 + y.^2 < disc_radius^2;
D = DictProfile(lines.ticks, func);


%% ===== Plot data ===== %%
figure();
bplines   = ScanLines(lines);
p0        = ProbeParams(0);
p0.angles = ProbeParam(angles);
bplines.params = p0;
bpimage = bplines.back_project();
subplot(131); lines.plot_lines();   title('Lines');
subplot(132); bpimage.draw_image(); title('Back project image');
subplot(133); lines.plot_psf();     title('Point-spread function');

pause();

%% ====== Signal Reconstruction ===== %%
% Generate algorithm 'IPalm' and solve the simulated problem.
p = deepcopy(lines.params);
p.intensity = ProbeParam(ones(lines.nlines,1),ones(lines.nlines,1)*[0.8,1.2],...
              @(v)ones(lines.nmeasures,1)*v');
lda = 0.05*max(pos(D*back_project(lines)));
prb = CalibLasso(lines,D,p,lda);
% alg = ReweightIPalmSecmRealdata(prb,D,lines);
alg = IPalmSecmRealdata(prb,D,lines);
alg.set_maxiter(500);
figure; 
alg.solve(); 
