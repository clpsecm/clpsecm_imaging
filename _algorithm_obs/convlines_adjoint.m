function LX = convlines_adjoint(RY,LD,direction)
% [Description]: 
%   Adjoint operator of D*L[X]
DIRECTION_FORWARD = 1;
DIRECTION_BACKWARD = -1;

[nbins,nangles] = size(RY);
d = length(LD);
dl = floor((d-1)/2);
dr = ceil((d-1)/2);


CDf = convmtx(LD,nbins);
CDr = convmtx(flip(LD),nbins);
LX = zeros(size(RY));
for I = 1:nangles
    RYI = [zeros(dl,1) ; RY(:,I) ; zeros(dr,1)];
    if direction(I) == DIRECTION_FORWARD
        LX(:,I) = CDf'*RYI;
    elseif direction(I) == DIRECTION_BACKWARD
        LX(:,I) = CDr'*RYI;
    end
end




