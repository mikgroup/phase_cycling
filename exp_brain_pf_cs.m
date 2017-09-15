%% Partial Fourier
clc
clear
close all
setPath

%% Get Kspace
load brain
ksp = ksp/max(max(max(abs(ifft2c(ksp)))));
FOV = size(ksp);
[sx, sy, nc] = size(ksp);

%% Undersampling
frac = 7/16; % partial fourier fraction
ax = 1.8; % sub-sampling factor in x
ay = 1.8; % sub-sampling factor in y
ncalib = 24;

mask = vdPoisMex(sx, sy, sx, sy, ax, ay, ncalib, 0, 2);
mask(:, 1:round(FOV(2)*frac)) = 0;
mask = repmat( mask, [1,1,FOV(3:end)]);
y = ksp.*mask;

mask = y ~= 0;

figure, imshow3( log(abs(y)) );
titlef('Log Magnitude of Kspace');


%% Generate Maps

ksize = [6, 6];
[maps, weights] = ecalib(y, ncalib, ksize);
figure, imshow3(abs(maps .* repmat( weights,1,1,nc)));
titlef('Sensitivity Maps');

%% Create linear operators

C = Identity;

S = ESPIRiT(maps, weights);

F = p2DFT(mask,[sx, sy, nc]);

M = Identity;

P = Identity;


%% Estimate Low Resolution Phase from calibration region

win = repmat(hanning(ncalib) * hanning(ncalib)', [1, 1, nc]);
p_est = angle(S' * (F' * zpad(crop(y, ncalib, ncalib, nc) .* win, sx, sy, nc)));

figure, imshow3(p_est)

% Incorporate phase in ESPIRiT linear operator
mapsp = maps .* repmat(exp(1j*p_est), [1, 1, nc]);
Sp = ESPIRiT(mapsp, weights);


%% Create proximal operators

lambdam = 0.003;
lambdap = 0.005;

Pm = wave_thresh('db4', 3, lambdam);

Pp = wave_thresh('db4', 3, lambdap);

%% Get Fully-sampled Image

x = S' * (ifft2c(ksp));

figure, imshowf(abs(x), [0, 1.0])
figure, imshowf(angle(x) .* (abs(x) > 0.05), [-pi, pi])


%% Zero-filled recon

x_sub = S' * ( F' * y );

figure, imshowf(abs(x_sub), [0, 1.0])
figure, imshowf(abs(abs(x_sub) - abs(x)), [0, 0.1])
figure, imshowf(angle(x_sub) .* (abs(x) > 0.05), [-pi, pi])

%% Bydder et al. Homodyne

maxiter = 1000;
doplot = 0;
dohogwild = 1;

lambdar = 0.003;
lambdai = 0.003;

Pr = wave_thresh('db4', 3, lambdar);
Pi = wave_thresh('db4', 3, lambdai);
tic
x_im = imhomodyne(y, Sp, F, Pm, Pp, maxiter, dohogwild, doplot);
toc
x_im = x_im .* sqrt(weights);

figure, imshowf(abs(x_im), [0, 1.0])
figure, imshowf(abs(abs(x_im) - abs(x)), [0, 0.2])
figure, imshowf(angle(x_im .* exp(1j * p_est)) .* (abs(x) > 0.1), [-pi, pi])

disp(psnr(abs(x), abs(x_im)))

%% Get phase wraps
ncycles = 16;
[m0, p0, W] = pfinit(y, S, F, ncycles);


%% Phase regularized reconstruction without phase cycling

niter = 100;
ninneriter = 10;
doplot = 1;
dohogwild = 1;

[mn, pn] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, {}, niter, ninneriter, dohogwild, doplot);

mn = mn .* sqrt(weights);

figure, imshowf(abs(mn), [0, 1.0])
figure, imshowf(abs(mn - abs(x)), [0, 0.2])
figure, imshowf(pn .* (abs(x) > 0.05), [-pi, pi])

disp(psnr(abs(x), abs(mn)))

%% Proposed phase regularized reconstruction with phase cycling

niter = 100;
ninneriter = 10;
doplot = 1;
dohogwild = 1;

[m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, W, niter, ninneriter, dohogwild, doplot);

m = m .* sqrt(weights);

figure, imshowf(abs(m), [0, 1.0])
figure, imshowf(abs(abs(m) - abs(x)), [0, 0.2])
figure, imshowf(p .* (abs(x) > 0.05), [-pi, pi])

disp(psnr(abs(x), abs(m)))

%% Zhao et al. Separate magnitude and phase reconstruction
% Requires irt from Jeff Fessler.
% Please run setup.m in the toolbox first.

lambda_m = 0.3; % regularization parameter for magnitude
lambda_p = 0.3; % regularization parameter for phase (rg2/rg4)

im_mask = weights > 0.1;
maps = maps .* im_mask;
proxg_m = wave_thresh('db4', 3, lambda_m);

y = ksp(mask == 1);
samp = mask(:, :, 1) == 1; % change to logical for irt.

[mi, xi] = separate_mag_phase_recon(y, samp, maps, im_mask, proxg_m, lambda_p);
disp(psnr(abs(x) / max(abs(x(:))), abs(mi) / max(abs(mi(:)))))

