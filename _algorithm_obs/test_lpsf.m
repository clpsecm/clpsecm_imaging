%% ===== Test 1: Fit Single PSF ====== %%
% %-Load Data
% blobfilename = 'blob_072816';
% fileType = '.xlsx';
% blobLength = 1500;
% 
% 
% blobDataRaw = xlsread([blobfilename, fileType]);
% blobLine = imresize(blobDataRaw(:,2),[blobLength,1]);
% blobLine = blobLine - min(blobLine); %-Adjust to be non-negative
% 
% 
% d = 50;
% LD = imresize(blobLine,[d,1]);
% LD = LD/sum(LD);
% 
% 
% %-Fit Curve
% p_init = [1e-1,2,-2,round(d/2),5e-2];
% p = nlinfit(1:numel(LD),LD',@lpsf,p_init);
% 
% 
% %-Plot Result
% idx = 1:numel(LD);
% plot(idx, LD,               'b-', ...
%      idx, lpsf(p,idx),      'r-', ...
%      idx, lpsf(p_init,idx), 'k-');



%% ===== Test 2: Fit One Line (with Gradient Descent) ====== %%
%---Line data---%
DATA_NUM = 2;
fileType = '.xlsx';
niter = 150;

switch DATA_NUM
    case 1
        filename = 'secm_sample_0202_072516';
        lineNumber = 2;
        lineLength = 4000; 
        xpos = [1260,1631,1905,2122,2506,2665];
    case 2
        filename = 'secm_sample_0104_071316';
        lineNumber = 1;
        lineLength = 3000;
        xpos = [1170 1425 1600 1900];
end
nbins = 200;
np = 5;
nx = numel(xpos);
d = nbins;

dataRaw = xlsread([filename, fileType] ,lineNumber);
dataLine = imresize(dataRaw(:,2),[lineLength,1]);


RY = imresize(dataLine,[nbins,1]);
RY = RY/sum(RY);
xpos_hat = round(xpos/lineLength*nbins);



%---Generate Kernel----%
%-Parameter
P = [5e-2,1.5,-1.5,round(d/2),2e-2]'*ones(1,nx);
P(4,:) = xpos_hat+round(5*rand(1,nx));
P = P(:);

%-Kernel
D = zeros(d,nx);
for I = 1:nx
    D(:,I) = lpsf(P((I-1)*np+1:I*np),1:d); 
end



%---Gradient Descent---%
Jp = zeros(nbins,np*nx);
tp = repmat([.01,.01,.01,100,.01]',nx,1);
dp = .1*P(1:np);
iiter = 1;
fig = figure;
fig1 = get(fig,'Number');
while iiter <= niter
    %-Jocobian
    for I = 1:nx
        for K = 1:np
            eK = zeros(np,1);
            eK(K) = 1;
            pK = P((I-1)*np+1:I*np);
            Jp(:,(I-1)*np+K) = ...
                (lpsf(pK+dp.*eK,1:d) - lpsf(pK,1:d)) / dp(K); 
        end
    end
    
    %-Gradient Descent
    Dp = zeros(nbins,nx);
    for I = 1:nx
        Dp(:,I) = lpsf( P((I-1)*np+1:I*np),1:d );
    end
    Res = sum(Dp,2) - RY;
    P = P - tp.*(Jp'*Res);   
    
    %-Show Result
    disp(['Residual = ', num2str(norm(Res))]);   
    disp(['Posotion = ', num2str( P(4+[0:nx-1]*np)' )]);
    figure(fig1); drawnow;
    subplot(121); plot(RY,'k') 
    subplot(122); plot(sum(Dp,2));
    
    %-Update iterator
    iiter = iiter+1;
end
