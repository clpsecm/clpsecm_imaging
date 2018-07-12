function D = inverse_Abel(LD)
%% [Description] 
%    Given the line projection data, assuming the data is circular 
%    symmetric, try to find the original 2D profile
%
%  [Algorithm]
%    Inverse Abel transform: 
%    D(r) = -1/pi int_r^inf ( dLD(y)/dy * 1/sqrt(y^2-r^2) ) dy

n = length(LD);
y = 0:(n-1);
r = y;


dLD = ([LD(2:end); 0] - LD);
Dr = zeros(n,1);
for I = 1:n %-r
    for J = I+1:n %-y
        Dr(I) = Dr(I) + dLD(J)/(sqrt(y(J)^2-r(I)^2));
    end
end
Dr = [-1/(pi)*Dr;0];


D = zeros(2*n-1,2*n-1);
for I = 1:2*n-1
    for J = 1:2*n-1
        dist = sqrt((I-n)^2 + (J-n)^2);
        K = dist+1;
        Ku = floor(K+1);
        Kd = floor(K);
        if K <= n
            D(I,J) = ((Ku-K)*Dr(Kd) + (K-Kd)*Dr(Ku));
        else
            D(I,J) = 0;
        end
    end
end


