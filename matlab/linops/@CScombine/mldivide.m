function res = mldivide(a,b)
% res = mtimes(FT, x)
%

pinvA = pinv( a.A );

ns = length(a.ppm);
[sx,sy,nc,nm,ne] = size( b );
assert( ne == length(a.TE) );


b = reshape( b, sx*sy*nc*nm, ne );
b = permute( b, [2, 1] );
res = pinvA * b;
res = permute(res, [2, 1] );
res = reshape( res, sx, sy, nc, nm, ns );
