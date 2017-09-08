function res = mtimes(a,b)

Nx = a.dataSize(1);
Ny = a.dataSize(2);
Nz = a.dataSize(3);
NTE = numel(a.TE);

if a.adjoint
    bb = reshape(b,[],NTE);     % prod(Nx,Ny,Nz) x NTE
    bb = permute(bb,[2 1]);     % NTE x prod(Nx,Ny,Nz)


    %this block implements (-j*2*pi*tn).*(conj(A*\rho_hat) .* bb
    rho_hat = permute(reshape(a.im_res,prod([Nx,Ny,Nz]),[]),[2 1]); % Nspecies x prod(Nx,Ny,Nz)
    res1 = conj(a.A*rho_hat) .* bb;        % NTE x prod(Nx,Ny,Nz)
    
    res1 = res1 .* ((-1i*2*pi*a.TE(:))*ones(1,prod([Nx,Ny,Nz])));     % NTE x prod(Nx,Ny,Nz)
    res1 = sum(res1,1);                                               % 1 x prod(Nx,Ny,Nz)
    
    %this block implements ctranspose(A)*bb
    res2 = bb;                      % NTE x prod(Nx,Ny,Nz)
    res2 = ctranspose(a.A) * res2;  % Nspecies x prod(Nx,Ny,Nz)
    
    res = cat(1,res1,res2);                         % Nspecies+1 x prod(Nx,Ny,Nz)
    res = reshape(permute(res,[2 1]),Nx,Ny,Nz,[]);  % Nx x Ny x Nz x Nspecies + 1
    
else
    bb = reshape(b,[],3);   % prod(Nx,Ny,Nz) x 3 (3 update terms, NOT TEs)
        
    % this block implements (A * \delta\rho)
    delta_rho = bb(:,2:3);                               % prod(Nx,Ny,Nz) x 2 (delta_water, delta_fat)
    res1 = a.A*permute(delta_rho,[2 1]);                 % NTE x prod(Nx,Ny,Nz)
    res1 = reshape(permute(res1,[2 1]),Nx,Ny,Nz,NTE);    % Nx x Ny x Nz x NTE
    
    % this block implements (j*2*pi*tn*\deltafield_map).*(A * \rho_hat)
    rho_hat = permute(reshape(a.im_res,prod([Nx,Ny,Nz]),[]),[2 1]); % Nspecies x prod(Nx,Ny,Nz)
    res2 = a.A*rho_hat;                                             % NTE x prod(Nx,Ny,Nz)
    tmp = bb(:,1);                                                  % prod(Nx,Ny,Nz) x 1
    res2 = res2 .* (1i*2*pi*a.TE(:)*transpose(tmp));                % NTE x prod(Nx,Ny,Nz) 
    res2 = reshape(permute(res2,[2 1]),Nx,Ny,Nz,NTE);               % Nx x Ny x Nz x NTE
    
    res = res1 + res2;  
end
