function [maps, weights] = ecalib( ksp, ncalib, ksize)

eigThresh_k = 0.03; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space

% Generate ESPIRiT Maps (Takes 30 secs to 1 minute)
[sx,sy,Nc] = size(ksp);
calib = crop(ksp,[ncalib,ncalib,Nc]);

% Get maps with ESPIRiT
[kernel,S] = dat2Kernel(calib,ksize);
idx = find(S >= S(1)*eigThresh_k, 1, 'last' );
[M,W] = kernelEig(kernel(:,:,:,1:idx),[sx,sy]);
maps = M(:,:,:,end);

% Weight the eigenvectors with soft-sense eigen-values
weights = W(:,:,end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;