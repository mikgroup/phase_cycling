function  res = B(A,field_map,TE,wf_est)

% Author:  Samir Sharma
% Created: September 2009

res.adjoint = 0;
res.imSize = [size(field_map,1) size(field_map,2) size(field_map,3)];
res.dataSize = res.imSize; 
res.A = A;
res.field_map = field_map; 
res.TE = TE;
res.im_res = wf_est;

res = class(res,'B');

