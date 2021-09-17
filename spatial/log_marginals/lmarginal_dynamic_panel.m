function results = lmarginal_dynamic_panel(y,x,W,N,T,info)
% PURPOSE: Bayesian log-marginal posteriors for space-time static and dynamic panel models: 
% static SLX, y(t) = alpha*iota + X(t)*b1 + W*X(t)*b2 + e(t) 
% static SDM, y(t) = alpha*iota + rho*W*y(t) + X(t)*b1 + W*X(t)*b2 + e(t)
% dynamic non-spatial panel (DLM), y(t) = alpha*iota + phi*y(t-1) + X(t)*b1 + W*X(t)*b2 + e(t)  
% dynamic (restricted) SDMR, y(t) = alpha*iota + phi*y(t-1) + rho*W*y(t) -rho*phi*W*y(t-1) + X(t)*b1 + W*X(t)*b2 + e(t)
% dynamic (unrestricted) SDMU, y(t) = alpha*iota + phi*y(t-1) + rho*W*y(t) + theta*W*y(t-1) + X(t)*b1 + W*X(t)*b2 + e(t)
%          b1, b2 =  no prior
%          rho,phi,theta = uniform prior
% static SDM sets
%-------------------------------------------------------------
% USAGE: results = lmarginal_dynamic_panel(y,x,W,N,T,info)
% where: y = dependent variable vector (N*T x 1)
%        x = independent variables matrix, WITH NO INTERCEPT TERM 
%        W = N by N spatial weight matrix
%        N = # of cross-sectional units
%        T = # of time periods plus 1
%        y = (T-1)*N
%        info.iflag = 0 for conventional W-matrix
%        info.iflag = 1 for transformed W-matirx (using dmeanF())
%        info.lflag = 0 for full lndet computation (default = 1, fastest)
%                   = 1 for MC lndet approximation (fast for very large problems)
%        info.order = order to use with info.lflag = 1 option (default = 50)
%        info.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%        info.rmin  = (optional) minimum value of rho to use in search (default = -1) 
%        info.rmax  = (optional) maximum value of rho to use in search (default = +1)    
%-------------------------------------------------------------
% RETURNS:  a structure:
%          results.meth   = 'lmarginal_dynamic_panel'
%          results.nobs   = # of cross-sectional observations
%          results.ntime  = # of time periods
%          results.nvar   = # of variables in x-matrix
%          results.y      = y-vector from input (N*T x 1)
%          results.time   = time for log-marginal posterior calculation
%          results.lflag  = lflag value from input (or default value used)
%          results.iflag  = iflag value from input (or default value used)
%          results.rmin   = minimum value of rho used (default +1)
%          results.rmax   = maximum value of rho used (default -1)
%          results.lmarginal = log marginal likelihood (5x1) [slx sdm dlm sdmr sdmu]'
%          results.probs     = model probabilities (5x1) [slx sdm dlm sdmr sdmu]'
% --------------------------------------------------------------
% NOTES: - returns only the log-marginal posterior (and probs) for model comparison purposes
%          NO ESTIMATES returned
% - results.lmarginal can be used for model comparison 
%
% this function uses bounds on the parameter space for integration
% that are less than the full range in an effort to speed computation
% specifically: the grid used is: -maxp*maxr - 0.5:incr:-maxp*maxr + 0.5;
% where maxp, maxr are posterior modes from the restricted model
% where theta = -rho*phi
% this usually works since theta is close to -rho*phi in most applications
% so the grid is defined on the basis of this
% you are of course free to broaden the grid using say:
% -maxp*maxr - 1:incr:-maxp*maxr + 1;
% --------------------------------------------------------------

% written by:
% James P. LeSage, last updated 12/2013
% Dept of Finance & Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

timet = clock; % start the timer

% error checking on inputs
[nt junk] = size(y);

[n1,nx] = size(x);

results.nobs  = N;
results.ntime = T;
results.nvar  = nx;
results.y = y;   

% set defaults
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
      error('lmarginal_dynamic_panel: Wrong # of input arguments');
     end;

results.rmin = rmin;
results.rmax = rmax;

pmin = -0.999;
pmax = 0.999;
tmin = -rmin*pmin;
tmax = 0.999;

IN = speye(N);

Tm1 = T-1;
e = ones(T-1,1);
L = spdiags(e, -1, T-1, T-1);

e = ones(T-1,1);
M = spdiags(e, 0, T-1, T-1);

% =========================================================
Wsmall = W;
% create large W-matrix for use later
Wb = sparse(kron(speye(Tm1),Wsmall)); 
xsdm = [ones(N*Tm1,1) x Wb*x];

% ====================================================================
% evaluate log-marginal for SLX model over a grid of rho values
xo = xsdm;
xpx = xo'*xo;

lndetx = log(det(xpx));
dof = (N*Tm1 -2*nx -1)/2;

% no uniform prior for rho in this model
logC = gammaln(dof) - dof*log(2*pi)  -0.5*lndetx;

bo = (xpx)\(xo'*y);
eo = y - xo*bo;
epeo = eo'*eo;

logm_slx = -dof*log(epeo) + logC;

results.logm_slx = logm_slx;
% ====================================================================
% evaluate log-marginal for SDM model over a grid of rho values

D = (1 - 1/rmin); % from uniform prior on rho

logC = -log(D) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx;

Wy = sparse(Wb)*y;
bo = (xpx)\(xo'*y);
bd = (xpx)\(xo'*Wy);
eo = y - xo*bo;
ed = Wy - xo*bd;
epeo = eo'*eo;
eped = ed'*ed;
epeod = ed'*eo;

incr = 0.001;
xxp=rmin:incr:rmax;
xx = xxp';
ngrid = length(xx);
iotan = ones(ngrid,1);
dof = (N*Tm1-2*nx-1)/2;

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
maxr = results.maxr_sdm;
% madj = adj - logm_sdm_profile;
yy = exp(logm_sdm_profile -adj);

% trapezoid rule integration
isum = trapz(xx,yy);
isum = isum + adj;

logm_out = isum + logC; % we put back the scale adjustment here
% 
results.logm_sdm = logm_out;

% ===========================================================
% do non-spatial dynamic lag model
% dof = (N*Tm1 -2*nx -1)/2;

% no uniform prior for rho in this model
% but a uniform prior over phi
% D = (1/pmax - 1/pmin) = 2, when pmax = 1, pmin = -1
logC = -log(2) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx;

xt = [x kron(M,W)*x ones(N*(T-1),1)];
xpx = xt'*xt;

% logdetx = log(det(xpx));

% logC = gammaln(dof) - dof*log(2*pi) -0.5*logdetx;

F1 = kron(M,IN)*y;
F2 = kron(L,IN)*y;

bo = (xpx)\(xt'*F1);
bp = (xpx)\(xt'*F2);

          E1 = (F1 - xt*bo);
          E2 = (F2 - xt*bp);

Q = zeros(2,2);
Q(1,1) = E1'*E1;
Q(1,2) = E1'*E2;
Q(2,1) = Q(1,2);
Q(2,2) = E2'*E2;

incr = 0.001;
xxp=pmin:incr:pmax;
xx = xxp';
ngrid = length(xx);
iotan = ones(ngrid,1);
qout = zeros(ngrid,1);

            qout = iotan*Q(1,1) - 2*xx*Q(1,2) + (xx.*xx)*Q(2,2);

            [adj,adin1] =  max(-dof*log(qout));
            
results.maxp = xx(adin1);
maxp = results.maxp;

            yy = exp(-dof*log(qout) - adj);

            isum = trapz(xx,yy);

            isum = isum + adj;

 logm_out = isum + logC; % we put back the scale adjustment here

results.logm_dlm = logm_out;

% do space-time SDM model with theta = -rho*phi
% =================================================================
% dof = (N*Tm1 -2*nx -1)/2;

% a uniform prior for rho and phi in this model
logC = -log(D) - log(2) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx;

xt = [x kron(M,W)*x ones(N*(T-1),1)];
xpx = xt'*xt;


F1 = kron(M,IN)*y;
F2 = kron(L,IN)*y;
F3 = kron(M,W)*y;
F4 = kron(L,W)*y;


bo = (xpx)\(xt'*F1);
bp = (xpx)\(xt'*F2);
br = (xpx)\(xt'*F3);
bt = (xpx)\(xt'*F4);

          E1 = (F1 - xt*bo);
          E2 = (F2 - xt*bp);
          E3 = (F3 - xt*br);
          E4 = (F4 - xt*bt);
       

Q = zeros(4,4);
Q(1,1) = E1'*E1;
Q(1,2) = E1'*E2;
Q(1,3) = E1'*E3;
Q(1,4) = E1'*E4;
Q(2,1) = Q(1,2);
Q(3,1) = Q(1,3);
Q(4,1) = Q(1,4);
Q(2,2) = E2'*E2;
Q(2,3) = E2'*E3;
Q(2,4) = E2'*E4;
Q(3,2) = Q(2,3);
Q(4,2) = Q(2,4);
Q(3,3) = E3'*E3;
Q(3,4) = E3'*E4;
Q(4,3) = Q(3,4);
Q(4,4) = E4'*E4;

% ====================================================================
% integrate log-marginal for space-time SDM model over a grid of rho, phi with theta constrained = -rho*phi values
incr = 0.005;
tmp=pmin:incr:pmax;
pgrid = tmp';
tmp = rmin:incr:rmax;
rgrid= tmp';
ngridp = length(pgrid);
ngridr = length(rgrid);

yy = sdm_fun(pgrid,rgrid,Q,W,dof,T,lflag,mcorder,mciter);

yadj = max(max(yy));

for i=1:ngridr;
    for j=1:ngridp;
        if yy(i,j) == yadj;
            maxr = rgrid(i,1);
            maxp = pgrid(j,1);
        end;
    end;
end;

% % [maxr maxp]

z1 = trapz(rgrid,exp(yy - yadj));
z3 = trapz(pgrid,z1);

% subplot(2,1,1),
% plot(z1);
% subplot(2,1,2),
% plot(z2);
% pause;

% z3
logm_out = z3 + yadj + logC;


results.logm_stsdm = logm_out;

% do space-time SDM model with theta unrestricted
% =================================================================
% ====================================================================
% integrate log-marginal for space-time SDM model over a grid of rho, phi and theta = -rho*phi values
% dof = (N*Tm1 -2*nx)/2;
% logC = -log(D) + gammaln(dof) - dof*log(2*pi) -0.5*logdetx;

tmp = -maxp*maxr - 0.5:incr:-maxp*maxr + 0.5;
tgrid = tmp';
ngridt = length(tgrid);

% a uniform prior for rho and
% a uniform prior for phi
% a uniform prior for theta 
logC = -2*log(D) - log(2) + gammaln(dof) - dof*log(2*pi)  -0.5*lndetx;


zsum = 0;
zout = zeros(ngridt,1);

for ii=1:ngridt; % loop over theta values and do bivariate integration on phi, rho
    theta = tgrid(ii);
    
yy = sdm_fun2(pgrid,rgrid,theta,Q,W,dof,T,lflag,mcorder,mciter);

yadj = max(max(yy));

z1 = trapz(rgrid,exp(yy - yadj));
z3 = trapz(pgrid,z1) + yadj;

zout(ii,1) = z3;

end;

yadj = max(max(zout));
% integrate out theta
z3 = trapz(tgrid,exp(zout - yadj));

% plot(tgrid,exp(zout - yadj));
% xlabel('theta');
% ylabel('theta posterior');
% pause;


% subplot(2,1,1),
% plot(z1);
% subplot(2,1,2),
% plot(z2);
% pause;

% z3

logm_out = z3 + yadj + logC;


results.logm_stsdmu = logm_out;

time = etime(clock,timet);

results.time = time;

% ===================================================
% calculate posterior model probabilities
lmarginal = [results.logm_slx results.logm_sdm results.logm_dlm results.logm_stsdm results.logm_stsdmu];

adj = max(lmarginal);
madj = lmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;

results.probs = probs';
results.lmarginal = lmarginal';
 

end
        

    function y = sdm_fun(pgrid,rgrid,Q,W,dof,T,lflag,mcorder,mciter)
            
            ngridp = length(pgrid);
            ngridr = length(rgrid);
          
            tmp = [];
            if lflag == 0
            lndet = lndetfull(rgrid,W);
            elseif  lflag == 1
            lndet = lndetmc(mcorder,mciter,W,rgrid);
            end
            for i=1:ngridp;
                tmp = [tmp lndet];
            end
            lndet = tmp;
            
       
            iotan = ones(ngridr,1);
            qout = zeros(ngridr,ngridp);
            
            for i=1:ngridp;
                p = pgrid(i,1);
                
                    r = rgrid;
                    t = -p*r;
                                
            qout(:,i) = iotan*Q(1,1) - 2*p*iotan*Q(1,2) + (p*p)*iotan*Q(2,2) - 2*r*Q(1,3) + (r.*r)*Q(3,3) -2*t*Q(1,4) + (t.*t)*Q(4,4) ...
                + (p*r)*Q(2,3) + (r*p)*Q(3,2) + (p*t)*Q(2,4) + (t*p)*Q(4,2) + (r.*t)*Q(3,4) + (t.*r)*Q(4,3);
            end;
                        
            y = -dof*log(qout) + T*lndet;
            
        end
    

    function y = sdm_fun2(pgrid,rgrid,theta,Q,W,dof,T,lflag,mciter,mcorder)
        
            
            ngridp = length(pgrid);
            ngridr = length(rgrid);
         
            tmp = [];
            if lflag == 0
            lndet = lndetfull(rgrid,W);
            elseif  lflag == 1
            lndet = lndetmc(mcorder,mciter,W,rgrid);
            end

            for i=1:ngridp;
                tmp = [tmp  lndet];
            end
            lndet = tmp;
                                 
            iotan = ones(ngridr,1);
            qout = zeros(ngridr,ngridp);
            
            for i=1:ngridp;
                p = pgrid(i,1);
                
                    r = rgrid;
                    t = theta;
                                
            qout(:,i) = iotan*Q(1,1) - 2*p*iotan*Q(1,2) + (p*p)*iotan*Q(2,2) - 2*r*Q(1,3) + (r.*r)*Q(3,3) -2*t*iotan*Q(1,4) + (t*t)*iotan*Q(4,4) ...
                + (p*r)*Q(2,3) + (r*p)*Q(3,2) + (p*t)*iotan*Q(2,4) + (t*p)*iotan*Q(4,2) + (r*t)*Q(3,4) + (t*r)*Q(4,3);
            end;
            
            
            y = -dof*log(qout) + T*lndet;
            
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


