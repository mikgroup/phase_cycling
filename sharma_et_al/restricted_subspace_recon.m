function [m, p] = restricted_subspace_recon(y, F, S, C, M, P, Pm, A, m, p, TE, niter, ninneriter, stepsize, minwinsize, damp)

% Start iteration
alpham = 1.0 / lipschitz(M);

Nx = size(p, 1);
Ny = size(p, 2);
NTE = size(y, 5);

current_support = [1 1];
XFM = Bspline(current_support, [Nx, Ny]);

for it = 1:niter
    
    expPp = exp(1j * (P * p));
    
    for inner_it = 1:ninneriter
        
        inner_it
        
        Mm = M * m;
        
        % Update residual, mag, and phase
        r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
        m = Pm(m + alpham * (M' * (conj(expPp) .* r )), alpham);
    figure(33)
    imshow3(abs(m))
    drawnow
    end
    
    % Update field map
    field_map = reshape(p(:, :, :, :, :, end), [Nx, Ny, 1]);
    residual = reshape(S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))), [Nx, Ny, 1, NTE]);
    
    e = update_field_map(residual, A, field_map, TE, squeeze(m), XFM);
    
    delta_p = real(e(:,:,1,1));
    p(:, :, 1, 1, 1, end) = p(:, :, 1, 1, 1, end) + damp * delta_p;
    
    if current_support == 1
        current_support = [Nx Ny];
    else
        current_support = max(round(stepsize * current_support), minwinsize);
    end
    figure(35)
    imshow3(p)
    drawnow
    current_support
    XFM = Bspline(current_support, [Nx, Ny]);
end