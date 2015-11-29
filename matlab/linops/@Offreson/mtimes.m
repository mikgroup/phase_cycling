function res = mtimes(a,b)
% res = mtimes(FT, x)
% out = in(1) + TE * in(2) 


if a.adjoint
    
    nf = 2;
    [sx,sy,nc,nm,ne] = size( b );
    assert( ne == length(a.TE) );
    
    
    b = reshape( b, sx*sy*nc*nm, ne );
    b = permute( b, [2, 1] );
    res = (a.A)' * b;
    res = permute(res, [2, 1] );
    res = reshape( res, sx, sy, nc, nm, nf );

else
    
    ne = length(a.TE);
    [sx,sy,nc,nm,nf] = size( b );
    assert( nf == 2 );
    
    b = reshape( b, sx*sy*nc*nm, nf );
    b = permute( b, [2, 1] );
    res = a.A * b;
    res = permute(res, [2, 1] );
    res = reshape( res, sx, sy, nc, nm, ne );
    
end

    

