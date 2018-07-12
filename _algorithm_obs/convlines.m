function RY = convlines(LX,LD,dir)
%% 
% [Description]:
%    Convolution on lines for asymmetric D
DIRECTION_FORWARD = 1;
DIRECTION_BACKWARD = -1;


[~,nangles] = size(LX);
d = length(LD);
dl = floor((d-1)/2);
dr = ceil((d-1)/2);

RY = zeros(size(LX));
for I = 1:nangles
    if dir(I) == DIRECTION_FORWARD
        RYI = conv(LX(:,I),LD);
    elseif dir(I) == DIRECTION_BACKWARD
        RYI = conv(LX(:,I),flip(LD));
    end
    RY(:,I) = RYI(dl+1:end-dr);
end




