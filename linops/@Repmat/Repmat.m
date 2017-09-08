function  res = Repmat(r)


res.adjoint = 0;
res.r = r;

res = class(res,'Repmat');

