function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint
    res = conj(a.m) .* b;
else
    res = a.m .* b;
end



    
