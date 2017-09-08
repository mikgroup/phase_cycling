function [m0, p0, W] = csinit(x, C, M, ncycles)

W = cell(ncycles,1);
[sx, sy, nc, nm, ne, ns] = size(C' * x);
A = getA(C);
mm = A \ reshape(permute(x, [5, 1, 2, 3, 4]), [ne, sx*sy*nc*nm]);
mm = permute(reshape(mm, [ns, sx, sy, nc, nm]), [2, 3, 4, 5, 6, 1]);
% mm = M \ (C' * x);
pp = x ./ ((C * (M * mm))+eps);
uu = pp(:, :, :, :, 2) .* conj(pp(:, :, :, :, 1));
m0 = abs(mm);

for c = 1:ncycles
    p = zeros(sx, sy, nc, nm, 1, ns+1);
    p(:, :, :, :, :, 1:ns) = angle(mm * exp(1j * (c-1) * 2 * pi /  ncycles));
    p(:, :, :, :, :, end) = angle(uu  * exp(1j * (c-1) * 2 * pi /  ncycles));
    
    if c == 1
        p0 = p;
    end
    
    W{c} = p - p0;
end
