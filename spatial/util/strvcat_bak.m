function f = strvcat(varargin)
% PURPOSE: replace char(varargin) with strvcat function
%          for compatability of newer matlab versions
%---------------------------------------------------
% USAGE:     vnames  = strvcat('name1','variable2','log(gdp)');
% where:     'name1','variable2','log(gdp)' = a series of variable length strings
%---------------------------------------------------
% RETURNS:
%        vnames = a vertical concatenation of fixed width strings
% --------------------------------------------------
% SEE ALSO: MATLAB char()
%---------------------------------------------------

nargs = length(varargin);

vnames = [];
for i=1:nargs;
vnames = [vnames 
          char(varargin{i})];
end;