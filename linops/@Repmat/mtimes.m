function res = mtimes(a,b)
% res = mtimes(FT, x)
%

if a.adjoint
    res = b;
    for i = 1:length(a.r)
        if a.r(i) > 1
            res = sum(res, i);
        end
    end
else
    res = repmat(b, a.r);
    
end

    
