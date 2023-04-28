
% Clear previous work
close all; 
clear; 
clc; 
addpath(genpath('data')); 
addpath(genpath('given')); 
%%%%%%%%%%%%%%%%% Load experimental data %%%%%%%%%%%%%%%%%

load('data/params.mat'); 
x_axis = params.x_axis; % Pixel locations along width [mm]
z_axis = params.z_axis; % Pixel locations along depth [mm]
Fs = params.Fs; % Sampling rate of the experiment

% Load the power-Doppler images
load('data/pdi.mat'); 
% PDI = PDI;

% Load the binary stimulus vector
load('data/stim.mat');
% stim = stim;

Nz = size(PDI,1); % Number of pixels along depth
Nx = size(PDI,2); % Number of width pixels
Nt = size(PDI,3); % Number of timestamps
t_axis = 0:1/Fs:(Nt-1)/Fs; % Time-axis of the experiment

%%%%%%%%%%%%%%%%%% Get to know the data %%%%%%%%%%%%%%%%%%

t_label = 100;
figure; imagesc(x_axis,z_axis,PDI(:,:,t_label)); 
ylabel('Depth [mm]'); xlabel('Width [mm]');
title(['PDI at ' num2str(round(t_label/Fs)) 'st second']);

figure; plot(t_axis,squeeze(PDI(10,10,:))); 
hold on; s = plot(t_axis,stim*3*10^6); % the stimulus is
  % multiplied with a factor to visualize it at the same 
                   % y-axis scale as the fUS time-series
xlabel('Time (s)'); ylabel('Power Doppler Amplitude');
title('Time-Series of the Pixel at (x,z)=(10,10)'); 
legend(s,'Stimulus');

% Calculate the mean PDI
mean_PDI = mean(abs(PDI), 3);
mean_PDI = mean_PDI./(max(mean_PDI(:)));
% Display log of mean_PDI to enhance the contrast
figure; imagesc(x_axis,z_axis,log(mean_PDI));  
title('Mean PDI');
ylabel('Depth [mm]'); xlabel('Width [mm]');

%%%%%%%%%%%%%%%%%%%%% Pre-processing %%%%%%%%%%%%%%%%%%%%%

% Standardize the PDI along time
P = (PDI - mean(PDI, 3))./std(PDI,[],3); 

% Spatial Gaussian smoothing
ht = fspecial('gaussian',[4 4],2);
Pg = double(convn(P,ht,'same'));

% Temporal low pass filter at 0.3 Hz per pixel time-series
f1 = 0.3;
[b, a] = butter(5,f1/(Fs/2),'low');
PDIlinear = reshape(Pg,Nz*Nx,Nt);
Pgf = reshape(filtfilt(b,a,PDIlinear')',size(PDI));
PDI = Pgf;

%%%%%%%%%%%%%%% Show the correlation image %%%%%%%%%%%%%%%
%%
% reshape PDI to 2-D
PCC = {};
PCC_aver = [];
PDIlinear = reshape(PDI,Nz*Nx,Nt);

% calculate max PCC
for lag= 1:10
    stim_lag = circshift(stim,lag);
    PCC{lag} = corr(PDIlinear.',stim_lag, 'Type','Pearson');
    PCC_aver(lag) = mean(abs(PCC{lag}));
end

lagindex = find(PCC_aver==max(PCC_aver)); 
%Note: also saved for question2 of CPD
PCC_sele = abs(PCC{ lagindex });
PCC_sele(find( abs(PCC_sele)<0.3 ))=0;
pc_image = reshape(PCC_sele, [Nz Nx]);

% Two ways to visualize the correlation image -->
plot_version = 1;
display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
    'Significantly Correlated Regions',plot_version);
% plot_version = 2;
% display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
%     'Significantly Correlated Regions',plot_version);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CPD %%%%%%%%%%%%%%%%%%%%%%%%%%%

% % You may use hidden_cpd_als_3d.m
% % Include plots for all your claims (you can use display_brain_img.m to 
% % help with the visualization of spatial maps)

R1 = 15;
corb3 = []; % correlation of temporal signature
options.maxiter = 300; 
options.th_relerr = 0.6;
% load('btd_r10_l2');
[B1, B2, B3, c, output] = cpd_als_3d(PDI, R1, options);

cpd_save={};
cpd_save{1} = B1;
cpd_save{2} = B2;
cpd_save{3} = B3;

%%
% Spatial Maps of CP Decomposed Components

% Uncomment following lines if fixed initialization needed
saveB = load('cpd_save2');
B1 = saveB.cpd_save{1};
B2 = saveB.cpd_save{2};
B3 = saveB.cpd_save{3};
figure;
for i = 1:size(B1,2)
    sc1 = B1(:,i)*B2(:,i)';
    % PCC > 0.3 shows srong correlation
    corb3(1,i) = corr(B3(:,i), circshift(stim,lagindex), 'Type','Pearson');
    subplot(5,3,i);
    imagesc(sc1);
    title(['corr(', num2str(i),')=',num2str(corb3(1,i))]);
    colorbar;
end
sgtitle('Spatial Maps of CP Decomposed Components');

%%
% Plot the relative error per iteration
figure; 
semilogy(output.relerr)
grid on
xlabel('Iteration number')
ylabel('Relative error $\frac{\| T-T_{dec} \|_F}{\| T\|_F}$','interpreter','latex');
title('Relative error of the CP decomposition')
%%
% Correlation between temporal signature and the stimulus
disp(['The strongest correlation appears at ' num2str(find(corb3>0.3 | corb3<-0.3))])
% manual choose the component index, here is 3
Corr_comp = corr(B3(:,3),circshift(stim,lagindex), 'Type','Pearson');

% The detailed spatial maps compared with correlation
% indices are manually selected as:
select1 = 6; select2 = 9; select3 = 10;
figure;
subplot(2,2,1); imagesc(pc_image); title('(a) correlation image');colorbar
subplot(2,2,2); imagesc(B1(:,select1)*B2(:,select1)');
title({['(b) component ', num2str(select1)]; ['corr=' num2str(corb3(1,select1))]}); colorbar
subplot(2,2,3); imagesc(B1(:,select2)*B2(:,select2)');
title({['(c) component ', num2str(select2)]; ['corr=' num2str(corb3(1,select2))]}); colorbar
subplot(2,2,4); imagesc(B1(:,select3)*B2(:,select3)');
title({['(d) component ', num2str(select3)]; ['corr=' num2str(corb3(1,select3))]}); colorbar
sgtitle({'Spatial Maps'; '- comparison between correlation and CPD component'});


%%
%%%%%%%%%%%%%%%%%%%%%%%%%% BTD %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fill in btd_ll1_als_3d.m
% Include plots for all your claims (you can use display_brain_img.m to 
% help with the visualization of spatial maps)

% Uncomment following lines if fixed initialization needed
% saveInit = load('12-2satisfied');
% init1 = saveInit.savedinitial;

R2 = 12;
L = 3;
corC = [];
options.maxiter = 300; 
options.th_relerr = 0.6;

% may use saved initialization as:
% z = load('12-2satisfied');
% init1 = z.savedinitial;
% [A, B, C, const, output,inioutput] = btd_ll1_als_3d(PDI, R2,L,options,init1);

% % or randomly use initials:
[A, B, C, const, output,inioutput] = btd_ll1_als_3d(PDI, R2, L, options);
savedinitial = inioutput;
%%
stim_lag = circshift(stim,lagindex);
figure;

for i = 1:R2
%     sc = A(:,2*i-1:(i*2))*B(:,i*2-1:i*2)'; % L=2
    sc = A(:,3*i-2:(i*3))*B(:,i*3-2:i*3)'; % L=3
    subplot(4,3,i);
    imagesc(sc);
    title(['corr(',num2str(i),')=',num2str(corr(C(:,i), stim_lag, 'Type','Pearson'))]);
    colorbar;
end
sgtitle('Spatial Maps of BT Decomposed Components');

%%
% Comparison of CPD and BTD
figure;

subplot(1,2,1);
imagesc(B1(:,10)*B2(:,10)');
title('Spatial map of CPD');
colorbar;
subplot(1,2,2);
imagesc(A(:,5:6)*B(:,5:6)');
title('Spatial map of BTD');
colorbar;
sgtitle('Comparison between CPD and BTD');
%%
% Plot the relative error per iteration
figure; 
semilogy(output.relerr)
grid on
xlabel('Iteration number')
ylabel('Relative error $\frac{\| T-T_{dec} \|_F}{\| T\|_F}$','interpreter','latex');
title('Relative error of the BT Decomposition')

%%

for i = 1:R2
    % PCC > 0.3 shows srong correlation
    corC(1,i) = corr(C(:,i), circshift(stim,lagindex), 'Type','Pearson');
end
disp(['Strong correlation appear at ' num2str(find(corC>0.3 | corC<-0.3))]);
disp(['The strongest correlation appears at ' num2str(find(max(abs(corC))))]);

% The detailed spatial maps compared with correlation maps
% indices are manually selected as:
select1 = 3; select2 = 7; select3 = 9;
% select1 = 8; select2 = 10; select3 = 11; %%BTD_QR
figure;
subplot(2,2,1); imagesc(pc_image); title('(a) correlation image');colorbar
subplot(2,2,2); imagesc(A(:,2*select1-1:(select1*2))*B(:,select1*2-1:select1*2)');
title({['(b) component ', num2str(select1)]; ['corr=' num2str(corr(C(:,select1), stim_lag, 'Type','Pearson'))]});
colorbar
subplot(2,2,3); imagesc(A(:,2*select2-1:(select2*2))*B(:,select2*2-1:select2*2)');
title({['(c) component ', num2str(select2)]; ['corr=' num2str(corr(C(:,select2), stim_lag, 'Type','Pearson'))]});
colorbar
subplot(2,2,4); imagesc(A(:,2*select3-1:(select3*2))*B(:,select3*2-1:select3*2)');
title({['(d) component ', num2str(select3)]; ['corr=' num2str(corr(C(:,select3), stim_lag, 'Type','Pearson'))]});
colorbar
sgtitle({'Spatial Maps'; '- comparison between correlation and BTD component'});
