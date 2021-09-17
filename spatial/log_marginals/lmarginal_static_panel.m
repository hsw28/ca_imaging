function results = lmarginal_static_panel(y,xo,W,N,T,info)
% PURPOSE: Bayesian log-marginal posterior for static spatial panel models
%          user should eliminate fixed effects using differencing, de-meaning transformations
%          no priors on beta, sige
%          uniform prior on rho, lambda over eigenvalue bounds
%-------------------------------------------------------------
% USAGE: results = lmarginal_static_panel(y,x,W,N,T,info)
% where: y = dependent variable vector (N*T x 1)
%        x = independent variables matrix, WITHOUT INTERCEPT TERM 
%        W = N by N spatial weight matrix (for W*y and W*e)
%        N = # of cross-sectional units
%        T = # of time periods
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       info.rmin  = (optional) minimum value of rho to use in search (default = -1) 
%       info.rmax  = (optional) maximum value of rho to use in search (default = +1)    
%       info.iflag = 0 for conventional W-matrix
%       info.iflag = 1 for transformed W-matirx (using dmeanF())
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth   = 'lmargainal_static_panel'
%          results.nobs   = # of cross-sectional observations
%          results.ntime  = # of time periods
%          results.y      = N*T x 1 vector of y from input
%          results.nvar   = # of variables in x-matrix
%          results.rmin   = minimum value of rho used (default +1)
%          results.rmax   = maximum value of rho used (default -1)
%          results.lflag  = lflag value from input (or default value used)
%          results.iflag  = iflag value from input (or default value used)
%          results.lmarginal = a 3 x 1 column-vector with [log-marginal]
%          [ logm_slx logm_sdm logm_sdem ]'
%          results.probs  = a 3 x 1 column-vector with model probs
%          results.logm_slx
%          results.logm_sdm
%          results.logm_sdem
% --------------------------------------------------------------
% NOTES: - returns only the log-marginal posterior and probabilities for model comparison purposes
%          NO ESTIMATES returned
% - results.lmarginal can be used for model comparison 
% --------------------------------------------------------------

% Koop (2003, p 42): "When comparing models using posterior odds ratios, it is
% acceptable to use noninformative priors over parameters which are common
% to all the models. However, informative, proper priors should be used over all
% other parameters."

% This function uses a uniform prior on rho, lambda, but no priors on beta,sigma
% 
% 
% written by:
% James P. LeSage, last updated 2/2014
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com


timet = clock; % start the timer

     rmin = -0.9999;
     rmax = 0.9999;
     mcorder = 30;
     mciter = 50;
     iflag = 0;
     lflag = 1;

     if nargin == 6
      fields = fieldnames(info);
      nf = length(fields);
      if nf > 0
       for i=1:nf
        if strcmp(fields{i},'order')
         mcorder = info.order;
        elseif strcmp(fields{i},'iter')
         mciter = info.iter;
        elseif strcmp(fields{i},'lflag')
         lflag = info.lflag;
         results.lflag = lflag;
        elseif strcmp(fields{i},'iflag')
         iflag = info.iflag;
         results.iflag = iflag;
         if iflag ==  0 % we have conventional W-matrix
          % default to using eigenvalue bounds for rho
          [n,junk] = size(W);
          opt.tol = 1e-4; opt.disp = 0;
          lambda = eigs(sparse(W),speye(n),1,'SR',opt);
          rmin = real(1/lambda) + 0.0001;
          rmax = 0.9999;
         elseif iflag == 1 % we have unconventional W-matrix
          lambda = eig(full(W));
          evalues = real(lambda);
          rmin = 1/min(evalues);
          rmax = 1/max(evalues);
         end
        end;
       end;
      end;
     elseif nargin == 5
      % set defaults
      mcorder = 30;
      mciter = 50;
      lflag = 1;
      rmin = -0.9999;
      rmax = 0.9999;
      results.lflag = 1;
      results.eig = 0;
      results.iflag = 0;
      
     else
      error('lmarginal_static_panel: Wrong # of input arguments');
     end;

results.rmin = rmin;
results.rmax = rmax;


[nt,nx] = size(xo);

results.nobs  = N;
results.ntime = T;
results.nvar  = nx;
results.y = y;   


IN = speye(N);
IT = speye(T);

Wsmall = W;
% create large W-matrix for use later
W = sparse(kron(speye(T),Wsmall)); 


xsdm = [ones(nt,1) xo W*xo];


% ====================================================================
% evaluate log-marginal for SDM model over a grid of rho values

x = xsdm;
xpx = x'*x;

lndetx_sdm = log(det(xpx));
if iflag == 0
 dof = (N*T -1)/2;
elseif iflag == 1
 Tm1 = T - 1;
 dof = (N*Tm1 -1)/2;
end;

D = (1 - 1/rmin); % from uniform prior on rho

logC_sdm = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sdm;

% logC_sdm
Wy = sparse(W)*y;

bo = (xpx)\(x'*y);
bd = (xpx)\(x'*Wy);
eo = y - x*bo;
ed = Wy - x*bd;
epeo = eo'*eo;
eped = ed'*ed;
epeod = ed'*eo;

% calculate log-marginal for SLX model
logC_slx = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx_sdm;

logm_slx = -dof*log(epeo) + logC_slx;

results.logm_slx = logm_slx;

incr = 0.001;
xxp=rmin:incr:rmax;
xx = xxp';
ngrid = length(xx);
iotan = ones(ngrid,1);

if lflag == 0
 if iflag == 0
  logm_sdm_profile = -dof*log(epeo*iotan - 2*xx*epeod + (xx.*xx)*eped) + T*lndetfull(xx,Wsmall);
  [adj,mind] = max(logm_sdm_profile);
 elseif iflag == 1
  logm_sdm_profile = -dof*log(epeo*iotan - 2*xx*epeod + (xx.*xx)*eped) + Tm1*lndetfull(xx,Wsmall);
  [adj,mind] = max(logm_sdm_profile);
 end;
else
 if iflag == 0
  logm_sdm_profile = -dof*log(epeo*iotan - 2*xx*epeod + (xx.*xx)*eped) + T*lndetmc(mcorder,mciter,Wsmall,xx);
  [adj,mind] = max(logm_sdm_profile);
 elseif iflag == 1
  logm_sdm_profile = -dof*log(epeo*iotan - 2*xx*epeod + (xx.*xx)*eped) + Tm1*lndetmc(mcorder,mciter,Wsmall,xx);
  [adj,mind] = max(logm_sdm_profile);
 end;
end;

results.maxr_sdm = xx(mind);
yy = exp(logm_sdm_profile -adj);

% trapezoid rule integration
isum = trapz(xx,yy);
isum = isum + adj;

logm_out = isum + logC_sdm; % we put back the scale adjustment here
% 
results.logm_sdm = logm_out;

% ======================================================================
% SDEM model

logC_sdem = -log(D) + gammaln(dof) - dof*log(2*pi);

% 
% 
% % do vectorized calculations

Wx = sparse(W)*x;
Wy = sparse(W)*y;
xpWx = x'*Wx;
xpWpx = Wx'*x;
xpWpWx = Wx'*Wx;
xpy = x'*y;
xpWy = x'*Wy;
xpWpy = Wx'*y;
xpWpWy = Wx'*Wy;
ypy = y'*y;
ypWy = y'*Wy;
ypWpy = Wy'*y;
ypWpWy = Wy'*Wy;


Q1 = zeros(ngrid,1);
Q3 = zeros(ngrid,1);
% 
for iter=1:ngrid;
    rho = xx(iter,1);
    
 Axx = xpx - rho*xpWx - rho*xpWpx + rho*rho*xpWpWx;
 
 Q3(iter,1) = log(det(Axx));
 
 Axy = xpy - rho*xpWy - rho*xpWpy + rho*rho*xpWpWy;

 Ayy = ypy - rho*ypWy - rho*ypWpy + rho*rho*ypWpWy;

 b = Axx\Axy;

 Q1(iter,1) = Ayy - b'*Axx*b;

end;

if lflag == 0
 if iflag == 0
  logm_sdem_profile = -dof*log(Q1) + T*lndetfull(xx,Wsmall) - 0.5*Q3;
  [adj2,mind] = max(logm_sdem_profile);
 elseif iflag == 1
  logm_sdem_profile = -dof*log(Q1) + Tm1*lndetfull(xx,Wsmall) - 0.5*Q3;
  [adj2,mind] = max(logm_sdem_profile);
 end;
else
 if iflag == 0
  logm_sdem_profile = -dof*log(Q1) + T*lndetmc(mcorder,mciter,Wsmall,xx) - 0.5*Q3;
  [adj2,mind] = max(logm_sdem_profile);
 elseif iflag == 1
  logm_sdem_profile = -dof*log(Q1) + Tm1*lndetmc(mcorder,mciter,Wsmall,xx) - 0.5*Q3;
  [adj2,mind] = max(logm_sdem_profile);
 end;
end
 

results.maxr_sdem = xx(mind);
yy = exp(logm_sdem_profile -adj2);

% trapezoid rule integration
isum = trapz(xx,yy);
isum = isum + adj2;

logm_out = isum + logC_sdem; % we put back the scale adjustment here

results.logm_sdem = logm_out;


% ===========================================================

time = etime(clock,timet);

results.time = time;

% calculate posterior model probabilities
lmarginal = [results.logm_slx results.logm_sdm results.logm_sdem];

adj = max(lmarginal);
madj = lmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;

results.probs = probs';
results.lmarginal = lmarginal';
 

end



function out=lndetfull(rvec,W)
% PURPOSE: computes Pace and Barry's grid for log det(I-rho*W) using sparse matrices
% -----------------------------------------------------------------------
% USAGE: out = lndetfull(W,lmin,lmax)
% where:    
%             W     = symmetric spatial weight matrix (standardized)
%             lmin  = lower bound on rho
%             lmax  = upper bound on rho
% -----------------------------------------------------------------------
% RETURNS: out = a structure variable
%          out.lndet = a vector of log determinants for 0 < rho < 1
%          out.rho   = a vector of rho values associated with lndet values
% -----------------------------------------------------------------------
% NOTES: should use 1/lambda(max) to 1/lambda(min) for all possible rho values
% -----------------------------------------------------------------------
% References: % R. Kelley Pace and  Ronald Barry. 1997. ``Quick
% Computation of Spatial Autoregressive Estimators'', Geographical Analysis
% -----------------------------------------------------------------------
 
% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

spparms('tight'); 
[n junk] = size(W);
z = speye(n) - 0.1*sparse(W);
p = colamd(z);
niter = length(rvec);
dettmp = zeros(niter,2);
for i=1:niter;
    rho = rvec(i);
    z = speye(n) - rho*sparse(W);
    [l,u] = lu(z(:,p));
    dettmp(i,1) = rho;
    dettmp(i,2) = sum(log(abs(diag(u))));
end;

out = dettmp(:,2);
end


function out=lndetmc(order,iter,wsw,xx)
% PURPOSE: computes Barry and Pace MC approximation to log det(I-rho*W)
% -----------------------------------------------------------------------
% USAGE: out = lndetmc(order,iter,W,rmin,rmax)
% where:      order = # of moments u'(wsw^j)u/(u'u) to examine (default = 50)
%              iter = how many realizations are employed (default = 30)
%                 W = symmetric spatial weight matrix (standardized)  
%              grid = increment for lndet grid (default = 0.01)
% -----------------------------------------------------------------------
% RETURNS: out = a structure variable
%          out.lndet = a vector of log determinants for -1 < rho < 1
%          out.rho   = a vector of rho values associated with lndet values
%          out.up95  = an upper 95% confidence interval on the approximation
%          out.lo95  = a lower  95% confidence interval on the approximation
% -----------------------------------------------------------------------
% NOTES: only produces results for a grid of 0 < rho < 1 by default
%        where the grid ranges by 0.01 increments
% -----------------------------------------------------------------------
% References: Ronald Barry and R. Kelley Pace, "A Monte Carlo Estimator
% of the Log Determinant of Large Sparse Matrices", Linear Algebra and
% its Applications", Volume 289, Number 1-3, 1999, pp. 41-54.
% -----------------------------------------------------------------------
 
 
% Written by Kelley Pace, 6/23/97 
% (named fmcdetnormgen1.m in the spatial statistics toolbox )
% Documentation modified by J. LeSage 11/99

[n,n]=size(wsw);

% Exact moments from 1 to oexact
td=full([0;sum(sum(wsw.^2))/2]);
oexact=length(td);

o=order;
% Stochastic moments

mavmomi=zeros(o,iter);
for j=1:iter;
u=randn(n,1);
v=u;
utu=u'*u;
for i=1:o;
v=wsw*v;
mavmomi(i,j)=n*((u'*v)/(i*utu));
end;
end;

mavmomi(1:oexact,:)=td(:,ones(iter,1));

%averages across iterations
avmomi=mean(mavmomi')';

clear u,v;

%alpha matrix

alpha=xx;
valpha=vander(alpha);
valphaf=fliplr(valpha);
alomat=-valphaf(:,(2:(o+1)));

%Estimated ln|I-aD| using mixture of exact, stochastic moments
%exact from 1 to oexact, stochastic from (oexact+1) to o

lndetmat=alomat*avmomi;

out = lndetmat;

end

