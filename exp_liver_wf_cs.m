%% Chemical Shift
clc
clear
close all
setPath

%% Get Kspace
load liver_waterfat

[sx,sy,nc,nm,ne] = size(ksp);
ns = length(ppm);
FOV = size(ksp);
ksp = ksp/max(vec(abs(ifft2c(ksp))));


%% Subsample
ax = 2; % sub-sampling factor in x
ay = 2; % sub-sampling factor in y
ncalib = 24;
mask = vdPoisMex(sx, sy, sx, sy, ax, ay, ncalib, 1, 1.3);
mask = repmat( mask, [1,1,nc, nm, ne]);
% mask = ones(size(mask));
figure, imshowf(mask(:, :, 1, :, 1))

y = ksp .* mask;

%% Generate Maps using ESPIRiT

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


%% Get initial solution and phase wraps
x = S' * (F' * y);
ncycles = 8;
[m0, p0, W] = csinit(x, C, M, ncycles);

figure, imshow3(m0)
figure, imshow3(p0)

%% Create proximal operators
lambdam = 0.003;
lambdap = 0.05;

Pm = wave_thresh('db4', 3, lambdam);

Pp = wave_thresh('db6', 3, lambdap);

%% Proposed phase regularized reconstruction with phase cycling

niter = 1000;
doplot = 1;
dohogwild = 1;

[m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m0, p0, W, niter, dohogwild, doplot);

figure, imshow3(m)
figure, imshow3(p)

%% Sharma et al.
A = calculate_chemical_shift_encoding_matrix(FieldStrength, ppm, TE);
lambdam = 0.01;
niter = 20;
ninneriter = 10;
stepsize = 0.75;
minwinsize = 16;
damp = 0.01;

Pm = wave_thresh('db4', 3, lambdam);
p0 = p0 * 0;
[m, p] = restricted_subspace_recon(y, F, S, C, M, P, Pm, A, m0, p0, TE, niter, ninneriter, stepsize, minwinsize, damp);


%%
figure, imshow3f(m)
sos = sqrt(sum(abs(m).^2, 6));
immask = sos > max(sos(:)) * 0.1;

figure, imshow3f(p .* immask)
