%% Phase Regularized Test
clc
clear
close all
setPath

%% Set Reconstruction Parameters
nIter = 500;
do_plot = true;


%% Get Kspace
load acetone

[sx,sy,nc,nm,ne] = size(k);
ns = length(ppm);

mask = k ~= 0;
k = k/max(max(max(abs(ifft2c(k)))));


figure, imshow3( log(abs(k)), [] );
titlef('Log Magnitude of Kspace');

%% Create linear operators

% create ESPIRiT operator
ESP_op = Identity;

% create partial Fourier operator
PFT_op = p2DFT(mask,[sx, sy, nc, nm, ne]);

% create water-fat combine operator
MAG_op = CScombine( TE, FieldStrength, ppm );

PHASE_op = Identity;
% PHASE_op = Offreson( TE, FieldStrength );

MAG_thresh = @(x, lambda) wave_thresh( x, lambda, 'db6' );

PHASE_thresh = @(x, lambda) wave_thresh( x, lambda, 'db8' );


%% Zero-filled recon
x_sub = ( PFT_op' * k );
figure,imshow3( cat( 2, abs(x_sub ) / max(abs(x_sub(:))) , abs( angle( x_sub ) ) / pi ),[] );


%% Separate magnitude and phase reconstruction 

MAGlambda = 0.001;
PHASElambda = 0.001;
MAGstep = 0.5;
PHASEstep = 1;


% Initialize mag and phase
[mag, phase] = sepMagPhaseRecon( k, ESP_op, PFT_op, ...
                                MAG_op, PHASE_op, ...
                                MAG_thresh, PHASE_thresh,...
                                MAGlambda, PHASElambda, ...
                                MAGstep, PHASEstep, ...
                                nIter, do_plot );

