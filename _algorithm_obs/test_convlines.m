

% <DRZ,D[RY,LD]> = <RY,D'[DRZ,LD]>
nbins = 100;
nlines = 10;
d = round(nbins/2);
dir = sign(randn(nlines,1));

RY = randn(nbins,nlines);
DRZ = randn(nbins,nlines);
LD = randn(d,1);


sum(sum(DRZ.*convlines(RY,LD,dir)))
sum(sum(RY.*convlines_adjoint(DRZ,LD,dir)))









