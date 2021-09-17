function stats=chainstats(chain)
%CHAINSTATS some statistics from the MCMC chain
% chainstats(chain,results)
%    chain    nsimu*npar MCMC chain
%    results  output from mcmcrun function

% $Revision: 1.4 $  $Date: 2009/08/13 15:47:35 $

mcerr = bmstd(chain)./sqrt(size(chain,1));

[z,p]  = geweke(chain);
tau    = iact(chain);
stats  = [mean(chain)',mcerr',tau', p'];



function [z,p]=geweke(chain,a,b)
%GEWEKE Geweke's MCMC convergence diagnostic
% [z,p] = geweke(chain,a,b)
% Test for equality of the means of the first a% (default 10%) and
% last b% (50%) of a Markov chain.
% See:
% Stephen P. Brooks and Gareth O. Roberts.
% Assessing convergence of Markov chain Monte Carlo algorithms.
% Statistics and Computing, 8:319--335, 1998.

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[nsimu,npar]=size(chain);

if nargin<3
  a = 0.1;
  b = 0.5;
end

na = floor(a*nsimu);
nb = nsimu-floor(b*nsimu)+1;

if (na+nb)/nsimu >= 1
  error('Error with na and nb');
end

m1 = mean(chain(1:na,:));
m2 = mean(chain(nb:end,:));

%%% Spectral estimates for variance
sa = spectrum0(chain(1:na,:));
sb = spectrum0(chain(nb:end,:));

z = (m1-m2)./(sqrt(sa/na+sb/(nsimu-nb+1)));
p = 2*(1-nordf(abs(z)));


function s=bmstd(x,b)
%BMSTD standard deviation calculated from batch means
% s = bmstd(x,b) - x matrix - b length of the batch
% bmstd(x) gives an estimate of the Monte Carlo std of the 
% estimates calculated from x

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:34 $

[n,p] = size(x);

if nargin<2
  b = max(10,fix(n/20));
end

inds = 1:b:(n+1);
nb = length(inds)-1;
if nb < 2
  error('too few batches');
end

y = zeros(nb,p);

for i=1:nb
  y(i,:)=mean(x(inds(i):inds(i+1)-1,:));
end

% calculate the estimated std of MC estimate
s = sqrt( sum((y - repmat(mean(x),nb,1)).^2)/(nb-1)*b );

function s=spectrum0(x)
%SPECTRUM0 Spectral density at frequency zero
% spectrum0(x) spectral density at zero for columns of x

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[m,n]= size(x);
s = zeros(1,n);
for i=1:n
  spec = spectrum(x(:,i),m);
  s(i) = spec(1);
end


function [y,f]=spectrum(x,nfft,nw)
%SPECTRUM Power spectral density using Hanning window
%  [y,f]=spectrum(x,nfft,nw) 

% See also: psd.m in Signal Processing Toolbox 

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:40 $

if nargin < 2 | isempty(nfft)
  nfft = min(length(x),256);
end
if nargin < 3 | isempty(nw)
  nw = fix(nfft/4);
end
noverlap = fix(nw/2);

% Hanning window
w = .5*(1 - cos(2*pi*(1:nw)'/(nw+1)));
% Daniel
%w = [0.5;ones(nw-2,1);0.5];
n = length(x);
if n < nw
    x(nw)=0;  n=nw;
end
x = x(:);

k = fix((n-noverlap)/(nw-noverlap)); % no of windows
index = 1:nw;
kmu = k*norm(w)^2; % Normalizing scale factor
y = zeros(nfft,1);
for i=1:k
% xw = w.*detrend(x(index),'linear');
  xw = w.*x(index);
  index = index + (nw - noverlap);
  Xx = abs(fft(xw,nfft)).^2;
  y = y + Xx;
end

y = y*(1/kmu); % normalize

n2 = floor(nfft/2);
y  = y(1:n2);
f  = 1./n*(0:(n2-1));


function y=nordf(x,mu,sigma2)
% NORDF the standard normal (Gaussian) cumulative distribution.
% NORPF(x,mu,sigma2) x quantile, mu mean, sigma2 variance

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.6 $  $Date: 2012/09/27 11:47:38 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

%y = 0.5*erf(-Inf,sqrt(2)*0.5*x);
y = 0.5+0.5*erf((x-mu)/sqrt(sigma2)/sqrt(2));


function [tau,m] = iact(dati)
%IACT estimates the integrated autocorrelation time
%   using Sokal's adaptive truncated periodogram estimator.

% Originally contributed by Antonietta Mira by name sokal.m

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if length(dati) == prod(size(dati))
  dati = dati(:);
end

[mx,nx] = size(dati);
tau = zeros(1,nx);
m   = zeros(1,nx);

x  = fft(dati);
xr = real(x);
xi = imag(x);
xr = xr.^2+xi.^2; %%%controllare questo
xr(1,:)=0;
xr=real(fft(xr));
var=xr(1,:)./length(dati)/(length(dati)-1);

for j = 1:nx
  if var(j) == 0
    continue
  end
  xr(:,j)=xr(:,j)./xr(1,j);
  sum=-1/3;
  for i=1:length(dati)
    sum=sum+xr(i,j)-1/6;
    if sum<0
      tau(j)=2*(sum+(i-1)/6);
      m(j)=i;
      break
    end
  end
end
