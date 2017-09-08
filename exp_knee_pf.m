%% Partial Fourier
clc
clear
close all
setPath

%% Get Kspace
load knee
ksp = ksp/max(max(max(abs(ifft2c(ksp)))));
FOV = size(ksp);
[sx, sy, nc] = size(ksp);

%% Undersampling
frac = 7/16;

r = 4;
ncalib = 20;
mask = zeros(sx, sy);
mask(:, 1:r:end) = 1;
mask(:, 1:round(FOV(2)*frac)) = 0;
mask(:, (sy-ncalib)/2:(sy+ncalib)/2) = 1;
figure, imshow3( mask);
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

%%

win = repmat(hanning(ncalib) * hanning(ncalib)', [1, 1, nc]);
p_est = angle(S' * (F' * zpad(crop(ksp, ncalib, ncalib, nc) .* win, sx, sy, nc)));

figure, imshow3(p_est)

mapsp = maps .* repmat(exp(1j*p_est), [1, 1, nc]);
Sp = ESPIRiT(mapsp, weights);


%% Create proximal operators

lambdam = 0.003;
lambdap = 0.01;

Pm = wave_thresh('db4', 3, lambdam);

Pp = wave_thresh('db4', 3, lambdap);

%% Get Fully-sampled Image

x = S' * (ifft2c(ksp));

figure, imshow(abs(x), [0, 1.0])
figure, imshow(abs(abs(x) - abs(x)), [0, 0.1])
figure, imshow(angle(x) .* (abs(x) > 0.1), [-pi, pi])


%% Zero-filled recon

x_sub = S' * ( F' * y );

figure, imshow(abs(x_sub), [0, 1.0])
figure, imshow(abs(abs(x_sub) - abs(x)), [0, 0.1])
figure, imshow(angle(x_sub) .* (abs(x) > 0.1), [-pi, pi])

%% Bydder

maxiter = 100;
doplot = 1;
dohogwild = 1;

x_im = imhomodyne(y, Sp, F, Pm, Pp, maxiter, dohogwild, doplot);
x_im = x_im .* sqrt(weights);

disp(psnr(abs(x), abs(x_im)))

figure, imshow(abs(x_im), [0, 1.0])
figure, imshow(abs(abs(x_im) - abs(x)), [0, 0.3])
figure, imshow(angle(x_im .* exp(1j * p_est)) .* (abs(x) > 0.1), [-pi, pi])


%% Get phase wraps
ncycles = 16;
[m0, p0, W] = pfinit(y, S, F, ncycles);


%% Proposed method without phase cycling

niter = 100;
doplot = 1;
dohogwild = 1;

[mn, pn] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, {}, niter, dohogwild, doplot);

mn = mn .* sqrt(weights);

disp(psnr(abs(x), abs(mn)))

%%
figure, imshow(abs(mn), [0, 1.0])
figure, imshow(abs(mn - abs(x)), [0, 0.3])
figure, imshow(pn .* (abs(x) > 0.1), [-pi, pi])

%% Proposed method with phase cycling

niter = 100;
doplot = 1;
dohogwild = 0;

[m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, W, niter, dohogwild, doplot);

m = m .* sqrt(weights);

disp(psnr(abs(x), abs(m)))

figure, imshow(abs(m), [0, 1.0])
figure, imshow(abs(abs(m) - abs(x)), [0, 0.3])
figure, imshow(p .* (abs(x) > 0.1), [-pi, pi])



