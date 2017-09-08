function [m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m, p, W, niter, dohogwild, doplot)

k = 0;
K = 10;
h = 1;

v = VideoWriter('iterations.mp4', 'MPEG-4');
v.FrameRate = 5
open(v)

% Start iteration
for it = 1:niter
    
    % Plot
    if (doplot)
        fig = figure(32);
        subplot(1,2,1),
        imshow3(abs(m))
        titlef(it);
        subplot(1,2,2),
        imshow3(real(p))
        titlef(it);
        writeVideo(v,getframe(fig));
        drawnow
    end
    
    Mm = M * m;
    expPp = exp(1j * (P * p));
    alpham = 1.0 / lipschitz(M) * h;
    alphap = 1.0 / lipschitz(P) / (max(abs(Mm(:)))^2 + eps) * h;
    
    % Get random phase wraps
    if isempty(W)
        w = 0;
    else
        t = randi(length(W));
        w = W{t};
    end
    
    % Update residual, mag, and phase
    r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
    m = Pm(m + alpham * real(M' * (conj(expPp) .* r )), alpham);
    p = Pp(p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r)))...
        + w, alphap) - w;
    
    if dohogwild
        k = k + 1;
        if k == K
            k = 0; 
            K = K * 2;
            h = h / 2;
        end
    end
    
end
close(v)