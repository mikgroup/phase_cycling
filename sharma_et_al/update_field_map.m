function error_est = update_field_map(residual,A,field_map,TE,wf_est,XFM)

%* BEGIN: FIELD MAP UPDATE *
[Nx,Ny,Nz,NTE] = size(residual);
Itnlim = 8;
FT = B(A,field_map,TE,wf_est);

residual = reshape(residual,[],NTE);    % prod(Nx,Ny,Nz) x NTE
residual_demod = residual .* exp(-1i*2*pi*conj(field_map(:))*TE); % prod(Nx,Ny,Nz) x NTE
residual_demod = reshape(residual_demod,Nx,Ny,Nz,NTE); % Nx x Ny x Nz x NTE

param = init;
param.FT = FT;
param.XFM = XFM;
param.data = residual_demod;
param.Itnlim = Itnlim;

res = zeros(Nx,Ny,Nz,3);

for n = 1:3
    res = fnlCg(res,param);
end
error_est = XFM'*res;
%* END: FIELD MAP UPDATE *

