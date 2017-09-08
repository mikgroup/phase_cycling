function x = imhomodyne(y, S, F, Pr, Pi, maxiter, dohogwild, do_plot)

x = S' * (F' * y);

k = 0;
K = 10;
h = 1;
for it = 1:maxiter
    r = (S' * (F' * (y - F * (S * x ))));
    x = x + h * r;
    x = Pr(real(x), h) + 1j * Pi(imag(x), h);
    
    
    if (do_plot)
        figure(32),
        subplot(1,2,1),
        imshow3( abs(x) )
        subplot(1,2,2),
        imshow3( abs(imag(x)) )
        titlef(it);
        drawnow
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