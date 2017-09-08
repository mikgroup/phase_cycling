function  res = lipschitz( op )

s = svd([eye(length(op.A)), op.A]);
res = s(1)^2;