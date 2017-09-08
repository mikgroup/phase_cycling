function res = mtimes(a,b)
% res = mtimes(FT, x)
% out = in(1) + TE * in(2) 


if a.adjoint
    
    [sx,sy,nc,nm,ne,ns] = size(b);
    assert( ne == length(a.TE) );
    
    res = zeros(sx, sy, nc, nm, 1, ns+1);
    res(:, :, :, :, :, 1:ns) = sum(b, 5);
    
    f = reshape(permute(b, [1, 2, 3, 4, 6, 5]), [sx*sy*nc*nm*ns, ne]);
    res(:, :, :, :, :, end) = sum(reshape(f * a.A, [sx, sy, nc, nm, 1, ns]), 6);

else
    
    ne = length(a.TE);
    ns = a.ns;
    [sx,sy,nc,nm,t,ns1] = size(b);
    assert(ns1 == ns+1);
    assert(t == 1);
    
    res = repmat(b(:, :, :, :, :, 1:ns), [1, 1, 1, 1, ne, 1]);
    f = repmat(b(:, :, :, :, :, end), [1, 1, 1, 1, 1, ns]);
    
    res = res + permute(reshape(f(:) * a.A', [sx, sy, nc, nm, ns, ne]), [1, 2, 3, 4, 6, 5]);
    
end

    

