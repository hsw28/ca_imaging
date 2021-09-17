% run cigarette demand model
% packs = f(tax,price,smokers,income)

clear all;

[a,b] = xlsread('pci.xls');
tn = 1:2:45;
income = a(2:end,tn);

incomem = income(:,2:end-5);

[a,b] = xlsread('states_borders.xls');

% creates a spatial weight matrix based on lengths of state borders in
% common with neighboring states
wraw = a(1:end,2:end);

n = 49;
wout = zeros(n,n);
for i=1:n;
    for j=1:n;
        wout(i,j) = wraw(i,j);
    end;
end;
 
for i=1:n;
    for j=1:n;
        wout(j,i) = wout(i,j);
    end;
end;
        
W = normw(wout); % normalizes the spatial weight matrix

load queenW.data;

Wc = normw(queenW); % 1st order contiguity weight matrix

T=17; % number of time periods
N=49; % number of regions

load state_tax_cig.data;
% data is stacked all states each year
%    col 1 = state tax revenue (in millions of dollars)
%    col 2 = tax per pack
%    col 3 = price per pack (including taxes)
%    col 4 = number of smokers

tax_rev = state_tax_cig(:,1);
tax_rev = tax_rev*1000000;
tax = state_tax_cig(:,2);
packs = tax_rev./tax;
price = state_tax_cig(:,3) - tax;
smokers = state_tax_cig(:,4);

y = log(packs); % dependent variable is packs of cigs sold legally

T=17; % number of time periods
N=49; % number of regions

ym = reshape(y,N,T);
taxm = reshape(tax,N,T);
pricem = reshape(price,N,T);
smokersm = reshape(smokers,N,T);

% truncate everything to deal with missing  values

y = ym(:,2:end);

taxt = taxm(:,2:end);
tax_lagt = taxm(:,1:end-1);

pricet = pricem(:,2:end);
price_lagt = pricem(:,1:end-1) - tax_lagt;

smokerst = smokersm(:,2:end);
smokers_lagt = smokersm(:,1:end-1);

incomet = incomem(:,2:end);
income_lagt = incomem(:,1:end-1);

yt = vec(y);

x1 = vec(log(taxt));

x2 = vec(log(pricet));

x3 = vec(log(smokerst));

x4 = vec(log(incomet));


xo = [x1 x2 x3 x4];



T=16; % number of time periods
N=49; % number of regions

% set ted=0 for model with spatial fixed effects without time dummies 
% set ted=1 for model with spatial and time period fixed effects
ted = 1;
[yf,xf,n,t,Wf]=demeanF(yt,xo,N,T,ted,W);

info.iflag = 1;

result1 = lmarginal_dynamic_panel(yf,xf,Wf,n,T,info); 


fprintf(1,'using W-miles weight matrix \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','dlm','sdmr','sdmu');
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal result1.probs];

mprint(out,in);

% using W-miles weight matrix 
% model log-marginal  model probs 
% slx       385.7177       0.0000 
% sdm       385.9187       0.0000 
% dlm       410.4034       0.4214 
% sdmr      410.1347       0.3221 
% sdmu      409.9068       0.2565 

% now use a different spatial weight matrix

[yf,xf,n,t,Wf]=demeanF(yt,xo,N,T,ted,Wc);

result2 = lmarginal_dynamic_panel(yf,xf,Wf,n,T,info); 

fprintf(1,'using W-contiguity weight matrix \n');
in.cnames = strvcat('log-marginal','model probs');
in.rnames = strvcat('model','slx','sdm','dlm','sdmr','sdmu');
in.width = 10000;
in.fmt = '%10.4f';
out = [result2.lmarginal result2.probs];

mprint(out,in);

% using W-contiguity weight matrix 
% model log-marginal  model probs 
% slx       387.2867       0.0000 
% sdm       388.1589       0.0000 
% dlm       412.5333       0.2998 
% sdmr      412.7585       0.3755 
% sdmu      412.6130       0.3247 

% compare miles versus contiguity matrix for SDEM models
tmpp = [result1.lmarginal 
        result2.lmarginal ];
        
adj = max(tmpp);
madj = tmpp - adj;
xx = exp(madj);
% compute posterior probabilities
psum = sum(xx);
probs1 = xx/psum;



in2.cnames = strvcat('W-miles','W-contiguity');
in2.rnames = strvcat('model','slx','sdm','dlm','sdmr','sdmu');
out = [probs1(1:5,1) probs1(6:10,1)];
        
mprint(out,in2);

% model      W-miles W-contiguity 
% slx         0.0000       0.0000 
% sdm         0.0000       0.0000 
% dlm         0.0329       0.2764 
% sdmr        0.0251       0.3462 
% sdmu        0.0200       0.2994 

% estimate sdmu model using W-contiguity matrix

% Use default prior settings

Wbig = kron(speye(t),Wf);

xnew  = [xf Wbig*xf];
k = size(xnew,2);

ndraw = 30000;
nomit = 20000;

results3 = sarstc_fe(yf,xnew,Wf,n,t,ndraw,nomit);

xnames = strvcat('tax','price','smokers','income');

prt_sarstc(results3,xnames);
