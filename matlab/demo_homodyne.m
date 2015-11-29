%% Phase Regularized Homodyne Test
clc
clear
close all
setPath

%% Set Reconstruction Parameters
nIter = 200;
do_plot = true;

MAGlambda = 0.005;
PHASElambda = 0.005;
MAGstep = 0.9;
PHASEstep = 0.9;

% ESPIRiT parameters
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % ESPIRiT kernel-window-size
eigThresh_k = 0.03; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space


%% Get Kspace
load knee
% load brain
k = k/max(max(max(abs(ifft2c(k)))));
FOV = size(k);


%% Undersampling
ax = 2;
ay = 2;
frac = 1/3;

mask = vdPoisMex( FOV(1), FOV(2), FOV(1), FOV(2), ax, ay, ncalib, 1, 2 );
mask(:, 1:round(FOV(2)*frac)) = 0;
mask = repmat( mask, [1,1,FOV(3:end)]);
k_sub = k.*mask;

mask = k_sub ~= 0;

figure, imshow3( log(abs(k_sub)), [] );
titlef('Log Magnitude of Kspace');


%% Generate Maps

[sx,sy,Nc] = size(k_sub);
calib = crop(k_sub,[ncalib,ncalib,Nc]);

% Get maps with ESPIRiT
[kernel,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_k));
[M,W] = kernelEig(kernel(:,:,:,1:idx),[sx,sy]);
maps = M(:,:,:,end);

% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;


figure, imshow3( abs(maps .* repmat( weights,1,1,Nc)), [] );
titlef('Sensitivity Maps');

%% Create linear operators

% create ESPIRiT operator
ESP_op = ESPIRiT(maps,weights);

% create partial Fourier operator
PFT_op = p2DFT(mask,[sx, sy, Nc]);

MAG_op = Identity;

PHASE_op = Identity;

MAG_thresh = @(x, lambda) db_wave_thresh( x, lambda );

PHASE_thresh = @(x, lambda) db_wave_thresh( x, lambda );


%% Get Fully-sampled Image

x = ESP_op' * ( ifft2c( k ) );
figure,imshow3( cat( 2, abs(x ) / max(abs(x(:))) , abs( angle( x ) ) / pi ) );
titlef('Fully-Sampled Image');

%% Zero-filled recon

x_sub = ESP_op' * ( PFT_op' * k_sub );
figure,imshow3( cat( 2, abs(x_sub ) / max(abs(x_sub(:))) , abs( angle( x_sub ) ) / pi ) );
titlef('Zero-filled Reconstruction');


%% Separate magnitude and phase reconstruction 

phase_init = angle(x_sub);

[mag, phase] = sepMagPhaseRecon( k_sub, phase_init, ESP_op, PFT_op, ...
                                MAG_op, PHASE_op, ...
                                MAG_thresh, PHASE_thresh,...
                                MAGlambda, PHASElambda, ...
                                MAGstep, PHASEstep, ...
                                nIter, do_plot );
