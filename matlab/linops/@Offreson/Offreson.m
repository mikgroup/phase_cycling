function  res = Offreson( TE, FieldStrength )


res.adjoint = 0;
res.TE = TE;

ne = length(TE);
nf = 2;

A = zeros(ne, nf);
GAMMA = 42.58; %in MHz/T

% df = GAMMA * FieldStrength * ppm(s);
for e = 1:ne
    A(e,1) =  1;
    A(e,2) =  TE(e);
end

s = svd(A);
A = A / s(1);

res.A = A;
res = class(res,'Offreson');

