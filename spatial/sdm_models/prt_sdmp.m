function prt_sdmp(results,vnames,fid)
% PURPOSE: Prints output using sdmp_g results structures
%---------------------------------------------------
% USAGE: prt_sdmp(results,vnames,fid)
% Where: results = a structure returned by sdm, sdm_g, sdm_gc, etc.
%        vnames  = an optional vector of variable names
%        fid     = optional file-id for printing results to a file
%                  (defaults to the MATLAB command window)
%--------------------------------------------------- 
%  NOTES: e.g. vnames = strvcat('y','const','x1','x2');
%         e.g. fid = fopen('ols.out','wr');
%  use prt_sdm(results,[],fid) to print to a file with no vnames               
% --------------------------------------------------
%  RETURNS: nothing, just prints the spatial regression results
% --------------------------------------------------
% SEE ALSO: prt, plt
%---------------------------------------------------   

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

if ~isstruct(results)
 error('prt_sdm requires structure argument');
elseif nargin == 1
 nflag = 0; fid = 1;
elseif nargin == 2
 fid = 1; nflag = 1;
elseif nargin == 3
 nflag = 0;
 [vsize junk] = size(vnames); % user may supply a blank argument
   if vsize > 0
   nflag = 1;          
   end;
else
 error('Wrong # of arguments to prt_sdm');
end;


nobs = results.nobs;
cflag = results.cflag;
p = results.p;
if cflag == 1
    nvars = p+1;
elseif cflag == 0
    nvars = p;
end;


% if (nflag == 1) % the user supplied variable names
% [tst_n nsize] = size(vnames);
%  if tst_n ~= nvars+1
%  fprintf(fid,'Wrong # of variable names in prt_sdm -- check vnames argument \n');
%  fprintf(fid,'will use generic variable names \n');
%  nflag = 0;
%  end
% end;
% 
% handling of vnames
% Vname = 'Variable';
% if nflag == 0 % no user-supplied vnames or an incorrect vnames argument
%     if cflag == 1 % a constant term
% 
%         Vname = strvcat(Vname,'constant');
%      for i=1:nvars-1
%         tmp = ['variable ',num2str(i)];
%         Vname = strvcat(Vname,tmp);
%      end;
%      for i=1:(nvars-1)
%         tmp = ['W*variable ',num2str(i)];
%         Vname = strvcat(Vname,tmp);
%      end;
%  
%     elseif cflag == 0 % no constant term
% 
%      for i=1:nvars
%         tmp = ['variable ',num2str(i)];
%         Vname = strvcat(Vname,tmp);
%      end;
%      for i=1:nvars
%         tmp = ['W*variable ',num2str(i)];
%         Vname = strvcat(Vname,tmp);
%      end;
%     end;
%  
%      
% % add spatial rho parameter name
%     Vname = strvcat(Vname,'rho');

    Vname = 'Variable';
    Vname = strvcat(Vname,vnames);
%      for i=1:nvars
%         Vname = strvcat(Vname,vnames(i+1,:));
%      end;
%      for i=1:nvars;
%         Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
%      end;
%     % add spatial rho parameter name
%         Vname = strvcat(Vname,'rho');
%      elseif cflag == 1 % a constant term
%      Vname = 'Variable';
%      for i=1:nvars
%         Vname = strvcat(Vname,vnames(i+1,:));
%      end;
%      for i=2:nvars;
%         Vname = strvcat(Vname,['W-' vnames(i+1,:)]);
%      end;
    % add spatial rho parameter name
        Vname = strvcat(Vname,'rho');
%     end; % end of cflag issue       
 
% end; % end of nflag issue



% find posterior means
    tmp1 = mean(results.bdraw);
    pout = mean(results.pdraw);
    bout = [tmp1'
        pout];
    tmp1 = std(results.bdraw);
    tmp2 = std(results.pdraw);
    bstd = [tmp1'
        tmp2];  


if strcmp(results.tflag,'tstat')
 tstat = bout./bstd;
 [junk nk] = size(results.bdraw);
 % find t-stat marginal probabilities
 tout = tdis_prb(tstat,results.nobs);
 results.tstat = bout./bstd; % trick for printing below
else % find plevels
 draws = [results.bdraw results.pdraw];
 [junk nk] = size(draws);
 for i=1:nk;
 if bout(i,1) > 0
 cnt = find(draws(:,i) > 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 else
 cnt = find(draws(:,i) < 0);
 tout(i,1) = 1 - (length(cnt)/(results.ndraw-results.nomit));
 end; % end of if - else
 end; % end of for loop
end; 



fprintf(fid,'\n');
fprintf(fid,'Bayesian Spatial Durbin Probit model\n');
if (nflag == 1)
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
end;
fprintf(fid,'# 0, 1 y-values    = %6d,%6d \n',results.zip,nobs-results.zip);
% fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',results.nobs,nk);
fprintf(fid,'ndraws,nomit       = %6d,%6d \n',results.ndraw,results.nomit);
fprintf(fid,'total time in secs = %9.4f   \n',results.time);
fprintf(fid,'time for sampling  = %9.4f \n',results.time3);

if results.lflag == 0
fprintf(fid,'No lndet approximation used \n');
end;
% put in information regarding Pace and Barry approximations
if results.lflag == 1
fprintf(fid,'Pace and Barry, 1999 MC lndet approximation used \n');
fprintf(fid,'order for MC appr  = %6d  \n',results.order);
fprintf(fid,'iter  for MC appr  = %6d  \n',results.iter);
end;
if results.lflag == 2
fprintf(fid,'Pace and Barry, 1998 spline lndet approximation used \n');
end;

fprintf(fid,'min and max rho= %9.4f,%9.4f \n',results.rmin,results.rmax);
fprintf(fid,'***************************************************************\n');
 if strcmp(results.tflag,'tstat')
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
      
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
 else % use p-levels for Bayesian results
tmp = [bout bstd tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Std Deviation'; pstring = 'p-level';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = Vname;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in);
end;


function bounds = cr_interval(adraw,hperc)
% PURPOSE: Computes an hperc-percent credible interval for a vector of MCMC draws
% --------------------------------------------------------------------
% Usage: bounds = cr_interval(draws,hperc);
% where draws = an ndraw by nvar matrix
%       hperc = 0 to 1 value for hperc percentage point
% --------------------------------------------------------------------
% RETURNS:
%         bounds = a 1 x 2 vector with 
%         bounds(1,1) = 1-hperc percentage point
%         bounds(1,2) = hperc percentage point
%          e.g. if hperc = 0.95
%          bounds(1,1) = 0.05 point for 1st vector in the matrix
%          bounds(1,2) = 0.95 point  for 1st vector in the matrix
%          bounds(2,1) = 0.05 point for 2nd vector in the matrix
%          bounds(2,2) = 0.05 point for 2nd vector in the matrix
%          ...
% --------------------------------------------------------------------

% Written by J.P. LeSage

% This function takes a vector of MCMC draws and calculates
% an hperc-percent credible interval
[ndraw,ncols]=size(adraw);
botperc=round((0.50-hperc/2)*ndraw);
topperc=round((0.50+hperc/2)*ndraw);
bounds = zeros(ncols,2);
for i=1:ncols;
temp = sort(adraw(:,i),1);
bounds(i,:) =[temp(topperc,1) temp(botperc,1)];
end;



