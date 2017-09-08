function res = init()
% param = init()
%
% function returns a structure with the entries that are needed for the reconstruction.
% The user MUST supply operators afterwords!
%
% See:
%	demo_SheppLoganTV demo_angio_simulation demo_Brain_2D.m
% 
% (c) Michael Lustig 2007

res.FT = [];                % measurement operator
res.XFM = [];               % sparse transform operator
res.data = [];              % measurements

res.Itnlim = 20;            % default number of iterations
res.gradToll = 1e-30;       % step size tollerance stopping criterea (not used)

% line search parameters
res.lineSearchItnlim = 150;
res.lineSearchAlpha = 0.01;
res.lineSearchBeta = 0.6;
res.lineSearchT0 = 1;       % step size to start with




