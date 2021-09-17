function [yf,xf,n,t,Wf]=demeanF(y,x,N,T,ted,W)
% set ted=0 for model with spatial fixed effects without time dummies 
% set ted=1 for model with spatial and time period fixed effects
n=N;
t=T;
nt=N*T;
if ted==0    % we do not have time dummy effects in the DGP
    Jt=speye(t)-1/t*ones(t,1)*ones(1,t);
    [Ftt junk]=eig(Jt);
    F=Ftt(:,2:t);
    y=reshape(y,n,t);
    y=y*F;
    y=reshape(y,n*(t-1),1);
    if isempty(x) == 0
        [junk,kx]=size(x);
        x=reshape(x,n,t,kx);
        xtemp=zeros(n,t-1,kx);
        for i=1:kx
            xtemp(:,:,i)=x(:,:,i)*F;
        end
        nt=nt-n;
        x=reshape(xtemp,nt,kx);
    else
        x=[];
    end
    t=t-1;
    Wf=W;
else % we have time dummy effects included in the DGP
    Jt=speye(t)-1/t*ones(t,1)*ones(1,t);
    [Ftt junk]=eig(Jt);
    F=Ftt(:,2:t);
    y=reshape(y,n,t);
    y=y*F;
    if isempty(x) == 0
        [junk,kx]=size(x);x=reshape(x,n,t,kx);xtemp=zeros(n,t-1,kx);
        for i=1:kx
            xtemp(:,:,i)=x(:,:,i)*F;
        end
        nt=nt-n;
        x=reshape(xtemp,nt,kx);
    else
        x=[];
    end
    t=t-1;
    
    Jn=speye(n)-1/n*ones(n,1)*ones(1,n);
    [Fnn junk]=eig(Jn);
    F=Fnn(:,2:n);
    y=F'*y;
    y=reshape(y,t*(n-1),1);   
    if isempty(x) == 0
        [junk,kx]=size(x);x=reshape(x,n,t,kx);xtemp=zeros(n-1,t,kx);
        for i=1:kx
            xtemp(:,:,i)=F'*x(:,:,i);
        end
        nt=nt-t;
        x=reshape(xtemp,nt,kx);
    else
        x=[];
    end
    n=n-1;
    Wf=F'*W*F;
end
yf=y;
xf=x;