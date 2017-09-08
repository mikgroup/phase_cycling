function [m0, p0, W] = pfinit(y, S, F, ncycles)

W = cell(ncycles,1);
x = S' * (F' * y);
m0 = abs(x);

for c = 1:ncycles
    p = angle(x * exp(1j * (c-1) * 2 * pi /  ncycles));
    
    if c == 1
        p0 = p;
    end
    
    W{c} = p - p0;
end
