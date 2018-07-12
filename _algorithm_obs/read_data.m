%% read_data.m
% [Description]: 
%  Read data from xls file, generate matlab data for recovery 
%  algorithm


%% ===== CONFIG ===== %%
mm = 1e-3;

DISPLAY_LINES = 1;
DISPLAY_ADJOINT = 1;
MODEL_TYPE = 'sym'; %'sym','asym'
IS_RAWDATA = 1;
LINE_INTEGRAL_METHOD = 'Fourier'; %-'Fourier','Direct'
DATA_NUMBER = 1;
IS_DLENGTH_USERDEF = 0;

pixelLength_mm = 30;
L = 0.35;   %-Scaling for lambda, lambda = L*DLt[RY]

if IS_DLENGTH_USERDEF == 1
    d = [21,21];
end

%% ===== Parameter Definition ===== %%

%-Global Variable
DIRECTION_FORWARD = 1;
DIRECTION_BACKWARD = -1;
fileType = '.xlsx';
if IS_RAWDATA
    dataTypeName = '_raw';
else
    dataTypeName = '';
end

%-Projected Blob Data File
blobFileName = 'blob_072816';
blobLength_mm = 1500;


%-Read Data File
switch DATA_NUMBER
    case 7 %-Data No.07 (Strange raw data: wrong offset)
    fileName   = ['secm_sample_0202_072516', dataTypeName];
    angles_deg = [0,30,60,90,150,180,240];
    offset_mm  = zeros(7,1);
    lineIdx   = [1:5];          
    lineLength_mm = 4000;  
    XsuppOpt   = [41.812, 41.67; 41.168, 41.391; 41.688, 40.604; 
                  41.591, 41.949; 42.765, 41.621; 42.17, 41.163];
    
    case 6 %-Data No.06 (Strange raw data: deg30 wrong, reversed line)
    fileName   = ['secm_sample_0104_071316', dataTypeName];
    angles_deg = [0, 30, 90, 120, 150, 180];
    offset_mm  = 1500 + [-400, 0, -600, -1400, -2000, -2800];
    lineIdx   = [1:5];
    lineLength_mm = 3000;
    XsuppOpt   = [24.498, 74.735; 25.036, 75.231;
                  25.207, 74.515; 25.303, 75.00];
    
    case 5 %-Data No.05 (Wrong Offset)
    fileName   = 'secm_sample_0202_071116';
    angles_deg = [0,30,60,90,120,150,180,240];
    offset_mm  = [0,0,500,1000,1000,1200];
    lineIdx   = [1:6];
    lineLength_mm = 4000; 
    XsuppOpt   = [41.812, 41.67; 41.168, 41.391; 41.688, 40.604; 
                  41.591, 41.949; 42.765, 41.621; 42.17, 41.163];
    
    case 4 %-Data No.04
    fileName   = ['secm_sample_0302_062116', dataTypeName];
    angles_deg = [0,15,30,45,60,75,90,105,120,135,165,180];
    offset_mm  = zeros(12,1);
    lineIdx   = [1:11];        
    lineLength_mm = 5000;

    case 3 %-Data No.03 (Heavily Skewed)
    fileName   = 'secm_sample_0302_060816';
    angles_deg = [0,15,30,45,60,90,105,120,150];
    offset_mm  = [1670,1300,1500,1300,1000,500,700,800,1370] - 2000;
    lineIdx   = [1:12];           
    lineLength_mm = 4000;

    case 2 %-Data No.02
    fileName   = ['secm_sample_0102_052616', dataTypeName];
    angles_deg = [0,30,60,90,120,150];
    offset_mm  = [1450, 800, 450, 300, 510, 900] - 1500;  
    lineIdx   = [1:6];
    lineLength_mm = 3000;
    XsuppOpt   = [24.538, 41.621; 25.23, 41.901; 25.412, 41.136];
    
    case 1 %-Data No.01
    fileName   = 'secm_sample_0302_012116';
    angles_deg = [0,15,30,45,60,75,90,105,120,135,150,165];
    endpoints  = [-1145,-1106,-992,-880,-723,-503,   0, 710,1373,1941,2377,2652;... 
                   0,    296, 573, 880,1251,1879,2345,2652,2377,1941,1373,710 ]';
    offset_mm  = sqrt(diag(endpoints*endpoints')) - 3500/2;
    lineIdx   = [1:12];
    lineLength_mm = 3500;
end


%% ===== Read data ===== %%
%---Select Subset of angles---%
nlines = length(lineIdx);
angles_deg = angles_deg(lineIdx);
offset_mm = offset_mm(lineIdx);


%---Line data---%
dataLines = zeros(lineLength_mm,nlines);
for I= 1:nlines
    dataRawI = xlsread([fileName, fileType] ,lineIdx(I));
    dataLines(:,I) = imresize(dataRawI(:,2),[lineLength_mm,1]);
end
clear dataRawI;


%-Flip lines, fix angle of slope fitting [-90°,90°]
angles_deg = mod(angles_deg,360);
lineDirection = zeros(nlines,1);
for I = 1:nlines
    if angles_deg(I) <= 90
        lineDirection(I) = DIRECTION_FORWARD;
    elseif angles_deg(I) > 90 && angles_deg(I) <= 270
        lineDirection(I) = DIRECTION_BACKWARD;
        angles_deg(I) = angles_deg(I) - 180;
        dataLines(:,I) = dataLines(end:-1:1,I);
        offset_mm(I) = -offset_mm(I);      
    elseif angles_deg(I) > 270 && angles_deg(I) <= 360
        lineDirection(I) = DIRECTION_FORWARD;
        angles_deg(I) = angles_deg(I) - 360;
    end 
end


%-[OPTIONAL] Read blob file
blobDataRaw = xlsread([blobFileName, fileType]);
blobLine = imresize(blobDataRaw(:,2),[blobLength_mm,1]);
blobLine = blobLine - min(blobLine); %-Adjust to be non-negative
clear blobDataRaw


%-Change unit for offset/slope angle
angles = angles_deg*pi/180;
offset = offset_mm/pixelLength_mm;


%-Define initial parameters for calibration
params_init.Amplitude = ones(1,nlines);
params_init.Rotation = angles_deg;
params_init.Translation = zeros(1,nlines); % [W.A] for symmetric
% params_init.Translation = offset;
params_init.Diameter = 100; 
params_init.Sigma = 80;
params_init.p = [0.1,1,-1];


%% ==== Process data ==== %%
%---Setup Dimenison of lines---%
nbins = floor(lineLength_mm/pixelLength_mm+1/2);

Ybd = get_domain_line_integral_adjoint(angles,offset,nbins);
Cy = [ (Ybd(1)+Ybd(2))/2 ; (Ybd(3)+Ybd(4))/2];

n = round([Ybd(2)-Ybd(1) , Ybd(4)-Ybd(3)]); % [ny,nx], dimension of Y


%---Downsampling---%
%-Reduce nbins by dropping samples, calibrate the size,
RY = zeros(nbins,nlines);
for I = 1 : nbins
    RY(I,:) = dataLines(floor(pixelLength_mm*(I-1/2)),:)*pixelLength_mm;
end
RY = RY - ones(nbins,1)*min(RY);
RY = RY./(ones(nbins,1)*sum(RY));


%---Generate Kernel Function---%
%-Profile D [TODO: Should Change to Function Handle]
if strcmp(MODEL_TYPE,'sym')
    %-Calculate diameter
    diam    = params_init.Diameter;
    sigma   = params_init.Sigma;
    gaudiam = ceil(4*sigma*2);
    if IS_DLENGTH_USERDEF == 0
        d = ceil((diam+gaudiam)/pixelLength_mm)*ones(1,2); %-dimension of D
    end
    
    %-Generate gaussian, disk
    dsk = fspecial('disk',diam/2);
    gau  = fspecial('gaussian',[gaudiam,gaudiam],sigma);

    %-Generate D
    D = imresize(conv2(dsk,gau),d);
    D = D/norm(D,'fro');

    %-Penalty variable lambda
    LtRY = fourier_line_integral_adjoint(RY,angles,n,offset,Cy);
    LDtRY = conv2_adjoint(LtRY,D);
    lambda = L*max(max(abs(LDtRY)));

    
elseif strcmp(MODEL_TYPE,'asym')
    %-Dimension of D
    d = floor(blobLength_mm/pixelLength_mm+1) * ones(1,2);
    didx = [1:d(1)]';
   
    %-Generate L[D]
    LDexp = imresize(blobLine,[d(1),1]);
    LDexp = LDexp/sum(LDexp);

    %-Generate D via Abel Inversion
    [~,maxBlobIdx] = max(blobLine);
    D = inverse_Abel(blobLine(maxBlobIdx:-1:1));
    D = imresize(D,d);
    D = D/norm(D,'fro');

    %-Find Fitting Parameters for assigned function class @lpsf from LD
    LD = @lpsf;
    p = nlinfit(didx,LDexp,LD,params_init.p);
    LDp = LD(p,didx);
    params_init.p = p;

    %-Penalty variable lambda
    DtRY = convlines_adjoint(RY,LDp,lineDirection);
    LDtRY = fourier_line_integral_adjoint(DtRY,angles,n,offset,Cy);
    lambda = L*max(max(abs(LDtRY)));
    
else 
    error('Wrong KERNEL_TYPE input');
end



%% ===== Plot data ===== %%
%---Adjoint Plot---%
if DISPLAY_ADJOINT   
    fig = figure;
    fig.WindowStyle = 'docked';
    LtRY = fourier_line_integral_adjoint(RY,angles,n,offset,Cy);
    imagesc(LtRY); axis equal; title('L^*[RY]')
end

%---Display Lines---%  
if DISPLAY_LINES
    fig = figure;
    fig.WindowStyle = 'docked';
    a = ceil(nlines/3);
    for I = 1:nlines;
        subplot(3,a,I); 
        plot(dataLines(:,I)); 
        title(['Raw Data: Scan ',num2str(angles_deg(I)),'^o']);
    end
    
    fig = figure;
    fig.WindowStyle = 'docked';
    for I = 1:nlines;
        subplot(3,a,I); 
        plot(RY(:,I)); 
        title(['Downsampled Data: Scan ',num2str(angles_deg(I)),'^o']);
    end
end


%-Dimension parameter
dim.n = n; 
dim.d = d;
dim.PixelLength = pixelLength_mm; 
dim.Offset = offset;
dim.Cy = Cy;
dim.Direction = lineDirection;


