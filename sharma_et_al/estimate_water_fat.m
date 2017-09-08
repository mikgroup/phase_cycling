function [wf_est,residual] = estimate_water_fat(image_comb,A,field_map,TE)

pinv_A = pinv(A);
[Nx,Ny,Nz,NTE] = size(image_comb);
image_comb = reshape(image_comb,[],NTE);            % prod(Nx,Ny,Nz) x NTE
image_comb = permute(image_comb,[2 1]);             % NTE x prod(Nx,Ny,Nz)
FM = exp(1i*2*pi*TE(:)*transpose(field_map(:)));    % NTE x prod(Nx,Ny,Nz)
rho_hat = pinv_A*(conj(FM).*image_comb);            % Nspecies x prod(Nx,Ny,Nz) 
residual = image_comb - FM .* (A*rho_hat);          % NTE x prod(Nx,Ny,Nz)
residual = reshape(permute(residual,[2 1]),Nx,Ny,Nz,NTE); % Nx x Ny x Nz x NTE
wf_est = reshape(permute(rho_hat,[2 1]),Nx,Ny,Nz,[]);     % Nx x Ny x Nz x Nspecies 


