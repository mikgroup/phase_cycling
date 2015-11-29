function res = mtimes(a,b)
% res = mtimes(FT, x)
%

if a.adjoint
    
    ns = length(a.ppm);
    [sx,sy,nc,nm,ne] = size( b );
    assert( ne == length(a.TE) );
    
    
    b = reshape( b, sx*sy*nc*nm, ne );
    b = permute( b, [2, 1] );
    res = (a.A)' * b;
    res = permute(res, [2, 1] );
    res = reshape( res, sx, sy, nc, nm, ns );

else
    
    ne = length(a.TE);
    [sx,sy,nc,nm,ns] = size( b );
    assert( ns == length(a.ppm) );
    
    b = reshape( b, sx*sy*nc*nm, ns );
    b = permute( b, [2, 1] );
    res = a.A * b;
    res = permute(res, [2, 1] );
    res = reshape( res, sx, sy, nc, nm, ne );
    
end

    
