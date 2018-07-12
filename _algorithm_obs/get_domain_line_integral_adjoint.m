function boundaryBox = get_domain_line_integral_adjoint(angles,offset,lineLength)
%%  
% [Description]: Given angle of slope for integrals and centers of line
%   integral, find the smallest squares containing it's adjoint.
%
% [Algorithm]: Adjoint of each line has domain as intersection of two
%   halfspaces, we aim to find the smallest square containing intersection
%   of all domain, with following steps
%   Obtain all hyperplains -> Find intersection -> Find boundary
% 
% [Input]:  angles(nangles):    Angles of slopes
%           offset(nangles):   Location of center on each lines from 
%                               middle in 0 degree (in pixel)
%           lineLength:         Scaler, length of line (in pixel)
% [Output]: boundaryBox(4,1):   [y1;y2;x1;x2] for each extreme point of square

%% ====== Parameter Definition ====== %%
n = length(angles); %-In this code we simplify nangles as n
boundaryBox = [-inf;inf;-inf;inf];
pL = -offset-lineLength/2; %-Left end points for each line
pR = -offset+lineLength/2;  %-Right end points for each line


%% ====== Find All Intersections ====== %%
% Intersect four lines by solving four linear equations:
%   x*cos(tI) + y*sin(tI) = pRI
%   x*cos(tI) + y*sin(tI) = pLI
%   x*cos(tJ) + y*sin(tJ) = pRJ
%   x*cos(tJ) + y*sin(tJ) = pLJ
for I = 1:(n-1)
    for J = (I+1):n
        bdIntxnPt = zeros(4,1); %-[y1;y2;x1;x2]
        invA = pinv([sin(angles(I)),cos(angles(I)); ...
                     sin(angles(J)),cos(angles(J))] );
        intxnPt = invA * [pL(I),pL(I),pR(I),pR(I) ; ...
                          pL(J),pR(J),pL(J),pR(J) ]; %-[py;px]
        
        bdIntxnPt(1) = min(intxnPt(1,:));
        bdIntxnPt(2) = max(intxnPt(1,:));
        bdIntxnPt(3) = min(intxnPt(2,:));
        bdIntxnPt(4) = max(intxnPt(2,:));
        
        boundaryBox(1) = max([boundaryBox(1),bdIntxnPt(1)]);
        boundaryBox(2) = min([boundaryBox(2),bdIntxnPt(2)]);
        boundaryBox(3) = max([boundaryBox(3),bdIntxnPt(3)]);
        boundaryBox(4) = min([boundaryBox(4),bdIntxnPt(4)]);
    end
end
