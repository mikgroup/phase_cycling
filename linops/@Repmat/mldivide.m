function res = mldivide(a,b)
% res = mtimes(FT, x)
%

res = b;
r = a.r;
for i = 1:length(r)
    if r(i) ~= 1
        res = sum(res, i) / r(i);
    end
end
