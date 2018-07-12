
%% ===== Test 1 ====== %%
% <RY,RX> = <Y,R'RX>
% n = [32,21];
% nbins = 80;
% nangles = 2;
% angles = pi*randn(1,nangles)-pi/2;
% % angles = [0,pi/2];
% Y = ones(n);
% RX = randn(nbins,nangles);
% 
% offset = [5,-15];
% centerY = [10;5];
% 
% 
% tic;
% RY = fourier_line_integral(Y,angles,nbins,offset,centerY);
% toc;
% tic;
% RtRX = fourier_line_integral_adjoint(RX,angles,n,offset,centerY);
% toc;
% 
% sum(sum(RY.*RX))
% sum(sum(Y.*RtRX))

%-This example gives supp(RY(:,1)) = {41:60}
%                    supp(RY(:,2)) = {22:51}                  


%% ===== Test 2 ===== %%
% row sum and column sum
% we can see some error induced 
% Y = ones(4,8);
% a = sum(Y)
% b = sum(Y')
% RY = fourier_line_integral(Y,[0,pi/2],8,[0;0],[0;0]);
% subplot(221); stem(a); title('sum(Y)');
% subplot(222); stem(b(end:-1:1)); title('sum(Y^t)');
% subplot(223); stem(RY(:,1)); title('RY, degree 0');
% subplot(224); stem(RY(:,2)); title('RY, degree 90');


%% ===== Test 3 ====== %%
% Point data
n = [100,100];
Y = zeros(n);
Y(20,40) = 1;
angles = [0;30;50;-40]*pi/180;
nangles = numel(angles);
nbins = 100;


RY = fourier_line_integral(Y,angles,nbins);
fig = figure;
fig.WindowStyle = 'docked';
subplot(411); stem(RY(:,1));
subplot(412); stem(RY(:,2));
subplot(413); stem(RY(:,3));
subplot(414); stem(RY(:,4));

fig = figure;
fig.WindowStyle = 'docked'; 
subplot(221); imagesc(fourier_line_integral_adjoint(RY(:,1),angles(1),n));
subplot(222); imagesc(fourier_line_integral_adjoint(RY(:,2),angles(2),n));
subplot(223); imagesc(fourier_line_integral_adjoint(RY(:,3),angles(3),n));
subplot(224); imagesc(fourier_line_integral_adjoint(RY(:,4),angles(4),n));



%% ===== Test 3 ==== %%
% Compare line_integral and fourier_line_integral
% Y = [zeros(40,51); ones(11,51)];
% nangles = 4;
% angles = pi*rand(nangles,1)-pi/2;
% % angles = [0,pi/2];
% RY1 = fourier_line_integral(Y,angles);
% RY2 = line_integral(Y,angles);
% for I = 1:nangles
%     subplot(3,nangles,I);         stem(RY1(:,I)); title(['Fourier ',num2str(I)]);
%     subplot(3,nangles,I+nangles); stem(RY2(:,I)); title(['Direct  ',num2str(I)]);
%     subplot(3,2,5); stem(sum(Y));
%     subplot(3,2,6); stem(sum(Y'));
% end







