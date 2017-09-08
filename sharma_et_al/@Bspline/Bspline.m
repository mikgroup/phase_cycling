function res = Bspline(current_support,full_res)

% Samir Sharma 
% June 2010

res.adjoint = 0;
[B,hx] = SDS_Create_BSpline(current_support,'cubic');
[res.basis] = SDS_BSpline_Basis(B,hx,full_res);
res = class(res,'Bspline');
