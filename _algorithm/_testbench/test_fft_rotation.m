% Generate square plot (odd pixel length)
N = 30;
M = ceil(N/2);
Y = zeros(N,N);
Y([M-6:M+6], [M-3:M+3]) = 1;

% Generate lines R from image Y
figure(); 
subplot(231); imagesc(Y);  set(gca,'YDir','normal'); title('0º'); set(gca,'YDir','normal');
subplot(232); imagesc(fft_rotate(Y,20));  set(gca,'YDir','normal'); title('20º');
subplot(233); imagesc(fft_rotate(Y,40));  set(gca,'YDir','normal'); title('40º');
subplot(234); imagesc(fft_rotate(Y,60));  set(gca,'YDir','normal'); title('60º');
subplot(235); imagesc(fft_rotate(Y,80));  set(gca,'YDir','normal'); title('80º');
subplot(236); imagesc(fft_rotate(Y,100));  set(gca,'YDir','normal'); title('100º');

% Generate squared plot (even pixel length)
N = 30;
M = N/2;
Z = zeros(N,N);
Z([M-5:M+6], [M-2:M+3]) = 1;

% Generate lines R from image Z
figure(); 
subplot(231); imagesc(Z);  set(gca,'YDir','normal'); title('0º');
subplot(232); imagesc(fft_rotate(Z,-30)); set(gca,'YDir','normal'); title('-30º');
subplot(233); imagesc(fft_rotate(Z,-60)); set(gca,'YDir','normal'); title('-60º');
subplot(234); imagesc(fft_rotate(Z,-90)); set(gca,'YDir','normal'); title('-90º');
subplot(235); imagesc(fft_rotate(Z,-120)); set(gca,'YDir','normal'); title('-120º');
subplot(236); imagesc(fft_rotate(Z,-150)); set(gca,'YDir','normal'); title('-150º');

% Test adjoints
% <Rt{Y1},Y2> = <R-t{Y2},Y1>
N = 51;
t = 60;
Y1 = randn(N,N);
Y2 = randn(N,N);
disp(['<Rt[Y1],Y2>  = ', num2str(sum(sum(fft_rotate(Y1, t).*Y2)))]);
disp(['<Y1,R-t[Y2]> = ', num2str(sum(sum(fft_rotate(Y2,-t).*Y1)))]);
