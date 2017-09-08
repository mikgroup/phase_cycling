function res = mtimes(a,b)

Nx = size(b,1);
Ny = size(b,2);
Nz = size(b,3);

res = zeros(Nx,Ny,Nz,3);
for aa = 1:numel(a.basis)
    support_size(aa) = sum(a.basis{aa}(1,:) > 0);
end

if a.adjoint
    for aa = 1:Nz
        if numel(support_size) == 1
            tmp = b(1:size(a.basis{1},1),1:size(a.basis{1},1),aa,1);
            tmp = transpose(a.basis{1}) * tmp;
            tmp = 1/support_size * transpose(a.basis{1}) * transpose(tmp);
        else
            tmp = b(1:size(a.basis{2},1),1:size(a.basis{1},1),aa,1);
            tmp = sqrt(1/support_size(2)) * transpose(a.basis{2}) * tmp;
            tmp = sqrt(1/support_size(1)) * transpose(a.basis{1}) * transpose(tmp);
        end
        res(:,:,aa,1) = tmp;
    end
else
    for aa = 1:Nz
        if numel(support_size) == 1
            tmp = a.basis{1} * b(:,:,aa,1);
            tmp = 1/support_size * a.basis{1} * transpose(tmp);
            res(1:size(a.basis{1},1),1:size(a.basis{1},1),aa,1) = tmp;
        else
            tmp = sqrt(1/support_size(1)) * a.basis{1} * b(:,:,aa,1);
            tmp = sqrt(1/support_size(2)) * a.basis{2} * transpose(tmp);
            res(1:size(a.basis{2},1),1:size(a.basis{1},1),aa,1) = tmp;
        end 
    end
end

res(:,:,:,2:3) = b(:,:,:,2:3);

end
