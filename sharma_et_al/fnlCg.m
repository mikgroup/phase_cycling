function [x] = fnlCg(x0,params)

% This is a modification of the function found in the SparseMRI
% toolbox - (c) Michael Lustig 2007
% Last Edited: Samir Sharma, November 2011

x = x0;

% line search parameters
maxlsiter = params.lineSearchItnlim;
gradToll = params.gradToll;
alpha = params.lineSearchAlpha; beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

g0 = wGradient(x,params);

dx = -g0;

while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx,FTXFMtdx] = preobjective(x,dx,params);
	f0 = objective(FTXFMtx,FTXFMtdx,0,params);
	t = t0;
    [f1, ERRobj, RMSerr]  =  objective(FTXFMtx,FTXFMtdx,t,params);
	
	lsiter = 0;
    
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 && (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx,FTXFMtdx,t,params);
    end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

 	% control the number of line searches by adapting the initial step search
 	if lsiter > 2
 		t0 = t0 * beta;
 	end 
 	if lsiter < 1
 		t0 = t0 / beta;
    end
    x = (x + t*dx);

 	%--------- uncomment for debug purposes ------------------------	
% 	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));
 	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	if (k > params.Itnlim) || (norm(dx(:)) < gradToll) 
		break;
	end
end
return;


function [FTXFMtx,FTXFMtdx] = preobjective(x,dx,params)
FTXFMtx = params.FT*(params.XFM'*x);
FTXFMtdx = params.FT*(params.XFM'*dx);


function [res,obj,RMS] = objective(FTXFMtx,FTXFMtdx,t,params)
obj = FTXFMtx + t*FTXFMtdx - params.data;
obj = obj(:)'*obj(:);
RMS = sqrt(obj/sum(abs(params.data(:))>0));
res = obj;


function grad = wGradient(x,params)
gradObj = gOBJ(x,params);
grad = gradObj;


function gradObj = gOBJ(x,params)
gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));
gradObj = 2*gradObj;
