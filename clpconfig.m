function clpcfg = clpconfig(ticks)
% clpcfg = CLPCONFIG(ticks) Configurate CLP
%
% [Userdef] User may change this function base on configuration of line
% probe. Currently there are four type of parameters for line probe
% {'angle', 'shifts', 'intensity', 'psf'}. 
%
% Each parameters  has its own operation on CLP probe:
%   'angles'    -  Determine the line scanning angle. 
%   'shifts'    -  Determine the shifts after line scan.
%   'intensity' -  Determine the the magnitude of per line, per dot 
%                  measurements after shifts.
%   'psf'       -  Determine the point spread function of line probe.
%
% Rules for configuration:
%   1. func(value) for has following formatting rules:
%       'angles'     -  func(value) is a vector(nlines,1)
%       'shifts'     -  func(value) is a vector(nlines,1)
%       'intensity'  -  func(value) is a matrix(nlines,nmeasures)
%       'psf'        -  func(value) is a vector(1,nmeasures)
%   2. 'bound' can be either n by 2 matrix or a generating function using
%      initial 'value'. The bound matrix consists of variable lower bound
%      (1st column) and upper bound (2nd column).
% * 3. If (param) is DEACTIVATED, set (param).value = NaN.
% * 4. If (param) is NOT A VARIABLE, set (param).bound =  @(v)v.
%
% See also PROBEPARAM, PROBEPARAMS
% disp(' ');
% disp('==== Load parameters from clpconfig =====')
% disp(' ');
if isrow(ticks); ticks = ticks'; end
nmeasures = numel(ticks);


%% ====== [Userdef] Probe Parameters ====== %%
% ----- Angles ----- %
clpcfg.angles.value = [0,45,70,90,115,135,180]';
clpcfg.angles.bound = @(v) v;
clpcfg.angles.func  = @(v) v;
  
% ----- Shifts ----- %
clpcfg.shifts.value = NaN;
clpcfg.shifts.bound = @(v) v;
clpcfg.shifts.func  = @(v) v;
 
% ---- Intensity ---- %
clpcfg.intensity.value = NaN;
clpcfg.intensity.bound = @(v) v; 
clpcfg.intensity.func  = @(v) ones(nmeasures,1)*v'; 

% ------- Psf ------- %
clpcfg.psf.value = [.6, .2, 5, 2, 1.5, 8]';
clpcfg.psf.bound = @(v) v;

t = ScanLines.get_psfticks(ticks);
fL = @(aL,cL) (t< 0).*min((-cL*t+1).^(-aL),realmax);
fR = @(aR,cR) (t>=0).*min(( cR*t+1).^(-aR),realmax);
fg = @(s) 1/sqrt(2*pi*s^2)*exp(-t.^2/(2*s^2));

clpcfg.psf.func = @(v) normalize(truncate(...
               conv( (1-v(1))*fL(v(3),v(4)) + v(1)*fR(v(5),v(6)), fg(v(2)) )...
               ,'mid') ,2);
% clpcfg.psf.func = @(v) circshift(normalize(truncate(...
%                conv( (1-v(1))*fL(v(3),v(4)) + v(1)*fR(v(5),v(6)), fg(v(2)) )...
%                ,'mid') ,2), v(7));

% ------------------- % 
           

 


