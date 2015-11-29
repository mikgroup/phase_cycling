function res = mtimes(a,x)
% This method applies the ESPIRiT Eigen-Vecs operator

maps = a.eigenVecs;

[sx,sy,nc,nv,ne] = size(maps);

if a.adjoint    
   %res = sum(conj(maps).*repmat(x,[1,1,1,nv]),3);    
    res = zeros(sx,sy,1,nv,ne);
    for n=1:nv
        res(:,:,n,:) = sum(conj(maps(:,:,:,n,:)).*x,3);
    end
else
   %res = sum(maps.*repmat(x,[1,1,nc,1]),4);   
   res = zeros(sx,sy,nc,1,ne);
   for n=1:nc
       res(:,:,n,:) =  sum(maps(:,:,n,:,:).*x,4);
   end
       
end

