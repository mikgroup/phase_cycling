function [m, p] = mprecon(y, F, S, C, M, P, Pm, Pp, m, p, W, niter, ninneriter, dohogwild, doplot)

k = 0;
K = 10;
h = 1;

% Start iteration
for it = 1:niter
    
    
    expPp = exp(1j * (P * p));
    alpham = 1.0 / lipschitz(M) * h;
    for itm = 1:ninneriter
    
        Mm = M * m;
        r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
        m = Pm(m + alpham * real(M' * (conj(expPp) .* r )), alpham);
        
        % Plot
        if (doplot)
            figure(32);
            subplot(1,2,1),
            imshow3(abs(m))
            titlef((it - 1) * ninneriter + itm);
            subplot(1,2,2),
            imshow3(real(p))
            titlef((it - 1) * ninneriter);
            drawnow
        end
    end
    
    alphap = 1.0 / lipschitz(P) / (max(abs(Mm(:)))^2 + eps) * h;
    Mm = M * m;
    for itp = 1:ninneriter
        
        % Get random phase wraps
        if isempty(W)
            w = 0;
        else
            t = randi(length(W));
            w = W{t};
        end
        
        expPp = exp(1j * (P * p));
        r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
        p = Pp(p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r)))...
            + w, alphap) - w;
        
        if (doplot)
            figure(32);
            subplot(1,2,1),
            imshow3(abs(m))
            titlef(it * ninneriter);
            subplot(1,2,2),
            imshow3(real(p))
            titlef((it - 1) * ninneriter + itp);
            drawnow
        end
    end
    
    if dohogwild
        k = k + 1;
        if k == K
            k = 0; 
            K = K * 2;
            h = h / 2;
        end
    end
    
end