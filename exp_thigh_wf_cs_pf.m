%% Chemical Shift
clc
clear
close all
setPath

%% Get Kspace
load thigh_waterfat

[sx,sy,nc,nm,ne] = size(ksp);
ns = length(ppm);
FOV = size(ksp);
ksp = ksp/max(vec(abs(ifft2c(ksp))));


%% Subsample
ax = 2; % sub-sampling factor in x
ay = 2; % sub-sampling factor in y
frac = 7/16; % partial fourier fraction
ncalib = 24;
mask = vdPoisMex(sx, sy, sx, sy, ax, ay, ncalib, 1, 1.3);
mask(:, 1:round(FOV(2)*frac)) = 0;
mask = repmat( mask, [1,1,nc, nm, ne]);
figure, imshowf(mask(:, :, 1, :, 1))

y = ksp .* mask;

%% Generate Maps

ksize = [6, 6];
[maps, weights] = ecalib(y(:, :, :, 1, 1), ncalib, ksize);
figure, imshow3(abs(maps .* repmat( weights,1,1,nc)));
titlef('Sensitivity Maps');

maps = repmat(maps, [1, 1, 1, 1, ne]);
weights = repmat(weights, [1, 1, 1, 1, ne]);

%% Create linear operators

C = ChemCombine(TE, FieldStrength, ppm);

S = ESPIRiT(maps, weights);

F = p2DFT(mask,[sx, sy, nc, nm, ne]);

M = Repmat([1, 1, 1, 1, ne, 1]);

P = Offres(TE, ns);


%% Get phase wraps
x = S' * (F' * y);
ncycles = 8;
[m0, p0, W] = csinit(x, C, M, ncycles);

figure, imshow3(m0)
figure, imshow3(p0)

%% Create proximal operators
lambdam = 0.003;
lambdap = 0.003;

Pm = wave_thresh('db4', 3, lambdam);

Pp = wave_thresh('db6', 3, lambdap);

%% Proposed phase regularized reconstruction with phase cycling

niter = 100;
ninneriter = 10;
doplot = 0;
dohogwild = 1;

tic
[m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, W, niter, ninneriter, dohogwild, doplot);
toc

figure, imshow3(m)
figure, imshow3(p)
