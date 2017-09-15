function [mi, xi] = separate_mag_phase_recon(y, samp, smap, mask, proxg_m, lamda_p)

%% Setup reconstruction parameters

nc = size(smap, 3);
Nd = [size(smap, 1), size(smap, 2)]; % image size
curv = prod(Nd); % spectral radius of A'*A

piter1 = 80; % # of iterations step 3 of initialization
piter2 = 80; % # of iterations for updating phase
nsubiter = [2 1]; % number of subiterations in each iteration

del = 0.005/1; % parameter for edge-preserving
preiter = 2; % number of iterations for intialization

%% Setup system objects
Q = Gdft('mask', mask, 'samp', samp, 'ifftshift', 1, 'fftshift', 1); % system fatrix

S = cell(nc,1);
F = cell(nc,1);
for i=1:nc
    s = smap(:,:,i);
    S{i} = Gdiag(s(mask), 'mask', mask);
    F{i} = Q;
end
S = block_fatrix(S, 'type', 'col');
F = block_fatrix(F, 'type', 'diag');
A = F * S;

%% Setup the wavlet transform object
% U = Godwt1(mask, 'level', 1);
% U = U';

%% The 1st step initialization by inverse FFT

img_i = A' * y / prod(Nd);
mi = abs(img_i);
xi = angle(img_i);

%% Regularized reconstruction: alternate updating for mt and xt

niter = piter1+piter2; % number of total iterations
C = Cdiffs([Nd(1),Nd(2)], 'mask', mask, 'order', 1);

for iter = 1:niter
    % update xi
    if iter > preiter
        if iter > piter1 % minimize the cost function
            xi = pcg_bls_exp_ep(A, C, y, mi, xi, lamda_p, del, ...
                nsubiter(1), 'mask', mask); % with regularizer 4
        else % optional 3rd step initialization
            xi = pcg_bls_ep(A, C, y, mi, xi, lamda_p, ...
                del, nsubiter, 'mask', mask);
        end
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % update mi or img_i by soft thresholding
    S = cell(nc,1);
    F = cell(nc,1);
    for i=1:nc
        s = smap(:,:,i);
        S{i} = Gdiag(s(mask) .* exp(1j * xi), 'mask', mask);
        F{i} = Q;
    end
    S = block_fatrix(S, 'type', 'col');
    F = block_fatrix(F, 'type', 'diag');
    Ax0 = F * S;
    
    if iter <= preiter
        % 2nd step initialization by conventional CS
        tmp = img_i + 1/curv * A' * (y - A * img_i);
        img_i = proxg_m(tmp, 1 / curv);
    else
        % CS update for magnitude image
        for subiter = 1:nsubiter(2)
            tmp = mi + 1/curv * real(Ax0'*(y-Ax0*mi));
            mi = proxg_m(tmp, 1 / curv);
            
        end
    end
    
    if iter == preiter
        mi = abs(img_i);
        xi = angle(img_i);
    end
    
    if mod(iter,10) == 1
        printf('iteration %d/%d', iter, niter)
    end
    
    figure(32);
    subplot(1,2,1),
    imshow3(abs(embed(mi,mask)))
    title(sprintf('Iteration: %d', iter));
    subplot(1,2,2),
    imshow3(real(embed(xi,mask)))
    title(sprintf('Iteration: %d', iter));
    drawnow
end

xi = embed(xi,mask);
mi = embed(mi,mask);

