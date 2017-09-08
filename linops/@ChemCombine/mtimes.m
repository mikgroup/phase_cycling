function res = mtimes(a,b)

if a.adjoint
    
    ns = length(a.ppm);
    [sx, sy, nc, nm, ne] = size(b);
    assert( ne == length(a.TE) );
    
    res = zeros(sx, sy, nc, nm, ne, ns);
    
    for s = 1:ns
        for e = 1:ne
            res(:, :, :, :, e, s) = res(:, :, :, :, e, s) + ...
                a.A(e, s)' * b(:, :, :, :, e);
        end
    end
else
    
    [sx, sy, nc, nm, ne, ns] = size(b);
    assert( ns == length(a.ppm) );
    
    res = zeros(sx, sy, nc, nm, ne);
    
    for s = 1:ns
        for e = 1:ne
            res(:, :, :, :, e) = res(:, :, :, :, e) + a.A(e, s) * b(:, :, :, :, e, s);
        end
    end
end

    
