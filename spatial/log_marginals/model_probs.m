function probs = model_probs(lmarginal)
% PURPOSE: computes and prints posterior model probabilities using log-marginals
% ---------------------------------------------------
%  USAGE: probs = model_probs(lmarginal)
%  where: log_marginal is an nmodel x 1 vector of scalars containing log
%  marginal likelihoods
%  returned by function, lmarginal_cross_section()
% e.g. result1 = lmarginal_cross_section(y,xo,W1);
%      lmarginal = result1.lmarginal(2);
%      result2 = lmarginal_cross_section(y,xo,W2);
%      lmarginal = [lmarginal
%                   result2.lmarginal(2)];
% model_probs(lmarginal);
% ---------------------------------------------------
%  RETURNS: probs = a vector of posterior model probabilities
% ---------------------------------------------------

% written by:
% James P. LeSage, 2/2014
% Dept of Finance & Economics
% Texas State University 
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com

% now scale using all of the vectors of log-marginals
% we must scale before exponentiating 
adj = max(max(lmarginal));
madj = lmarginal - adj;

xx = exp(madj);

% compute posterior probabilities
psum = sum(xx);
probs = xx/psum;

