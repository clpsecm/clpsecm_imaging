% Test <R,L[Y]> = <Lt[R],Y>

% Basinc setting
ticks = [-50:50];
nmeasures = numel(ticks);
nlines    = 4;

% Setup ProbeParams
params = ProbeParams(NaN);
angles = 180*rand(nlines,1);
shifts = 10*randn(nlines,1);
mask   = randn(nmeasures,nlines);
kernel = randn(nmeasures,1);
params.angles    = ProbeParam(angles,angles,@(v)v);
params.shifts    = ProbeParam(shifts,shifts,@(v)v);
params.intensity = ProbeParam(1, 1,  @(v) mask   );
params.psf       = ProbeParam(2, 2,  @(v) normpdf(ticks,0,v)' );
      
% Random imge Y
image = randn(nmeasures,nmeasures);
Y = SecmImage(ticks,image);

% Random lines R
currents = randn(nmeasures,nlines);
R = ScanLines(ticks,currents,params);

% The following should be equal
lhs = inprod(R.back_project(),       Y)
rhs = inprod(Y.line_project(params), R)