function  res = Offres(TE, ns)

% TE = TE / norm(TE);

res.adjoint = 0;

% Make TE(1) = 0 and TE(2) = 1 for simplicity
TE = TE - TE(1);
TE = TE / TE(2);
res.TE = TE(:);
res.ns = ns;

A = TE(:);

res.A = A;
res = class(res,'Offres');