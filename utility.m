% The functions provided here support the feasibility of the main code file.
% Note that the MFE toolbox is required to run this code "https://www.kevinsheppard.com/MFE_Toolbox".
% Y.ZHAO Apr/2018

%% The LambdaOmega.m compute the CUSUM statistics proposed by Aue et al (2009). Function required: Sk.m, vechdata.m, LRcov.m
function [lambda, omega, ind, mat]=LambdaOmega(data,kernel_type)
[n,~]=size(data);
[sk] = Sk(data);
[vech]=vechdata(data);
switch lower(kernel_type)
    case{'bartlett'}
        [Cov] = LRcov(vech','bartlett');
    case{'parzen'}
        [Cov] = LRcov(vech','parzen');
    case{'truncated'}
        [Cov] = LRcov(vech','truncated');
    case{'tukey'}
        [Cov] = LRcov(vech','tukey');
    case{'qs'}
        [Cov] = LRcov(vech','qs');
    case{'flat'}
        [Cov] = LRcov(vech','flat');
    otherwise
        error('The kernel type must be "Bartlett", or "Parzen", or "Tukey", or "qs", or "tf", or ''flat''.');
end
   mat=zeros(n,1);
   for i=1:n
       mat(i,1)=sk(:,i)'*inv(Cov)*sk(:,i);
   end
   [lambda,ind]=max(mat);
   omega=sum(mat)/n;
end


%% SK.m function compute the partial sum statistics Sk
function [Sk] = Sk(matrix)
[n,d]=size(matrix);
ybar=sum(matrix)/n;
newmat=NaN(n,d);
omg=d*(d+1)/2;
Sk=NaN(omg,n);
for i=1:n
    newmat(i,:)=matrix(i,:)-ybar;
end
vechmat=NaN(omg,n);
for i=1:n
    vechmat(:,i)=vech(newmat(i,:)'*newmat(i,:));
end
for i=1:n
    sk=(sum(vechmat(:,1:i),2)-(i/n)*sum(vechmat(:,:),2))/sqrt(n);
    Sk(:,i)=sk;
end
end

%% vechdata.m provides the vech operator, function required: vech.m
function [ vechdata ] = vechdata( data )
[n,d]=size(data);
vechdata=zeros(d*(d+1)/2,n);
for i=1:n
    vechdata(:,i)=vech(data(i,:)'*data(i,:));
end
end
% vech.m
function stackedData = vech(matrixData)
[k,l] = size(matrixData);
sel = tril(true(k));
stackedData = matrixData(sel);
end

%% LRcov.m compute the long run covariance matrix using different types of kernel functions
% Matlab function required: CUSUMtest.m, optbandw.m, scovnw.m
function [V] = LRcov(data,kernel_type)
[n,d]=size(data);
Ind=zeros(d,1);
StatM=zeros(d,1);
for i=1:d
    [Stat,~,ind]=CUSUMtest(data(:,i));
    StatM(i,1)=Stat;
    Ind(i,1)=ind;
end
dmdata=zeros(size(data));
for i=1:d
    dmdata(1:Ind(i,1),i)=data(1:Ind(i,1),i)-mean(data(1:Ind(i,1),i));
    dmdata((Ind(i,1)+1):n,i)=data((Ind(i,1)+1):n,i)-mean(data((Ind(i,1)+1):n,i));
end
switch lower(kernel_type)
    case{'bartlett'}
        bw = optbandw(dmdata(:,1),'bartlett');
        V=scovnw(dmdata,bw,1,'hacc_b');
    case{'parzen'}
        bw = optbandw(dmdata(:,1),'parzen');
        V=scovnw(dmdata,bw,1,'hacc_p');
    case{'truncated'}
        bw = optbandw(dmdata(:,1),'tf');
        V=scovnw(dmdata,bw,1,'hacc_t');
    case{'tukey'}
        bw = optbandw(dmdata(:,1),'tukey');
        V=scovnw(dmdata,bw,1,'hacc_th');
    case{'qs'}
        bw = optbandw(dmdata(:,1),'qs');
        V=scovnw(dmdata,bw,1,'hacc_qs');
    case{'flat'}
        bw = n^(0.25);
        V=scovnw(dmdata,bw,1,'hacc_f');
    otherwise
        error('The kernel type must be "Bartlett", or "Parzen", or "Tukey", or "qs", or "tf", or ''flat''. ');
end
end

%% CUSUMtest.m function [Stat, stat, ind] = CUSUMtest(series)
% one function performs non-parametric CUSUM test to detect change point.
function [Stat, stat, ind] = CUSUMtest(series)
n=length(series);
Tn=zeros(n,1);
for i=1:n
    tn=(sum(series(1:i,1))-(i/n)*sum(series))/sqrt(n);
    Tn(i,1)=tn;
end
Xkhat=zeros(n,1);
Xkstar=zeros(n,1);
Var=zeros(n,1);
for i=1:n
    xkhat=sum(series(1:i,1))/i;
    Xkhat(i,1)=xkhat;
    xkstar=sum(series((i+1):n,1))/(n-i);
    Xkstar(i,1)=xkstar;
end
for i=1:n-1
    var=(sum((series(1:i,1)-Xkhat(i,1)).^2)+sum((series((i+1):(n-1),1)-Xkstar(i,1)).^2))/n;
    Var(i,1)=var;
end
Var(n,1)=(sum((series-Xkhat(n,1)).^2))/n;
stat=zeros(n,1);
t=1:n;
t=t./n;
w=sqrt(t.*(1-t));
w=w';
for i=1:n
    stat(i,1)=1/sqrt(Var(i,1))*abs(Tn(i,1))./w(i,1);
end
[Stat,ind]=max(stat);
end

%% optbandw.m compute the optimal bandwidth (Newey and Wests's Method)
function bandwidth = optbandw(moments, kernel_type)
[T,q] = size(moments);
% START WITH ERROR CHECK
if nargin<1, error('Insufficient number of inputs.'), end;
if T==1, error('You need more than on observation. Check the size of the moments.'), end;
if nargin<2, kernel_type = 'bartlett', disp('The default method (Bartlett) is used'),end;

% The main procedure starts here
b = moments*ones(q,1);

switch lower(kernel_type)
    case{'bartlett'}
        nu = 1;
        n = fix((T^(1/9))^2);
        cgamma = 1.1447;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';      % NW bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       
        j = (1:n)';                            % NW bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));  
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1));% NW bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
    case{'parzen'}
        nu  = 2;
        n   = fix((T^(1/25))^4);
        cgamma=2.6614;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';            % parzen bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       % |
        j = (1:n)';                                  % |----> parzen bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));        % |
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1)); % parzen bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
    case{'tukey'}
        nu  = 2;
        n   = fix((T^(1/5))^4);
        cgamma=1.7462;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';            % tukey bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       % |
        j = (1:n)';                                  % |----> tukey bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));        % |
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1)); % tukey bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
    case{'qs'}
        nu  = 2;
        n   = fix((T^(1/25))^4);
        cgamma=1.3221;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';            % qs bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       % |
        j = (1:n)';                                  % |----> qs bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));        % |
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1)); % qs bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
    case{'tf'}
        nu  = 2;
        n   = fix((T^(1/25))^4);
        cgamma=1.3221;
        c = toeplitz(b,b(1:n+1));
        transfb = repmat(b,1,n+1);
        sigma = (1/T)*(sum(c.*transfb))';            % tf bandwidth selection, step 2
        s0 = sigma(1,1)+2*sum(sigma(2:end,1));       % |
        j = (1:n)';                                  % |----> tf bandwidth selection, step 3 
        snu = 2*sum((j.^nu).*sigma(2:end,1));        % |
        gammahat = cgamma*((snu/s0)^2)^(1/(2*nu+1)); % tf bandwidth selection, step 4
        bandwidth = fix(gammahat*(T^(1/(2*nu+1))));
        bandwidth = 0.2266*bandwidth;
        otherwise
        error('The kernel type must be "Bartlett", or "Parzen", or "Tukey", or "qs", or "tf".');
end

%% scovnw.m calculates the Long-run covariance estimation using Newey-West (Bartlett) weights 
function V=scovnw(data,nlag,demean,type)
T=size(data,1);
if nargin==1
    nlag=min(floor(1.2*T^(1/3)),T);
    demean=true;
elseif nargin==2
    demean=true;    
end    
if isempty(nlag)
    nlag=min(floor(1.2*T^(1/3)),T);
end
if isempty(demean)
    demean=true;
end
if ~ismember(demean,[0 1]) 
    error('DEMEAN must be either logical true or false.')
end
if floor(nlag)~=nlag || nlag<0 
    error('NLAG must be a non-negative integer.')
end
if ndims(data)>2
    error('DATA must be a T by K matrix of data.')
end
if demean
    data=data-repmat(mean(data),T,1);
end
% NW weights
if strcmp(type, 'hacc_b')==1
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=1-a;
% Truncated weights
elseif strcmp(type, 'hacc_t')==1 
w=ones(1,nlag+1);
% Parzen weights
elseif strcmp(type, 'hacc_p')==1 
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag+1
       if a(1,i)<=0.5
          w(1,i)=1-6*a(1,i)^2+6*abs(a(1,i))^3;
       else
          w(1,i+1)=2*(1-abs(a(1,i)))^3;
       end
   end
% Tukey-Hanning weights
elseif strcmp(type, 'hacc_th')==1 
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag+1
       w(1,i)=(1+cos(pi*a(1,i)))/2;
   end  
% Quadratic Spectral weights
elseif strcmp(type, 'hacc_qs')==1 
   a=[1/(nlag):1/(nlag):1];
   a=[0,a];
   w=ones(1,nlag+1);
   for i=1:nlag
       w(1,i+1)=(25/(12*pi^2*a(1,i+1)^2))*(sin(6*pi*a(1,i+1)/5)/(6*pi*a(1,i+1)/5)-cos(6*pi*a(:,i+1)/5));
   end
% Flat Top weights
elseif strcmp(type, 'hacc_flat')==1
    a=[1/(nlag):1/(nlag):1];
    a=[0,a];
    w=ones(1,nlag);
    for i=1:nlag
        if a(1,i)<=0.05
        w(1,i)=1;
        else
        w(1,i)=exp(-1*exp(-1/(a(1,i)-0.05)^2)/(a(1,i)-1)^2);
        end
    end
   w=[w,0];
end
V=data'*data/T;
for i=1:nlag
    Gammai=(data((i+1):T,:)'*data(1:T-i,:))/T;
    GplusGprime=Gammai+Gammai';
    V=V+w(i+1)*GplusGprime;
end
end

%% Constant Conditional Correlation data generating process
function [re]=concorrsim(cor,numob)
 wnumob=numob+0.2*numob;% warming-up parameter 0.2
 epsilon=randn(2,wnumob);
 R=[1,cor;cor,1];
 [Q,~] = chol(R);% Cholesky factorization on conditional correlation matrix at each time period to keep positive-definate property on forthcoming variance-covariance matrix.
 Uchol = Q';
 Pchol = Uchol*Uchol';
 for i=1:wnumob
 epsilon(:,i) = Uchol*epsilon(:,i);% Let generated innovation series has covariance structure as Uchol. 
 end
epsilon=transpose(epsilon);

% conditional variance processes
 re1=zeros(wnumob,1);
 re2=zeros(wnumob,1);
 h1=zeros(wnumob,1);
 h2=zeros(wnumob,1);
 h1(1,1)=1;
 h2(1,1)=1;
 for i=1:(wnumob-1);
    re1(i,1)=sqrt(h1(i,1))*epsilon(i,1);
    h1(i+1,1)=0.01+0.01*(re1(i,1)*re1(i,1))+0.94*h1(i,1);
 end
 re1(wnumob,1)=sqrt(h1((wnumob-1),1))*epsilon((wnumob-1),1);
 for i=1:(wnumob-1);
    re2(i,1)=sqrt(h2(i,1))*epsilon(i,2);
    h2(i+1,1)=0.01+0.01*(re2(i,1)*re2(i,1))+0.94*h2(i,1);
 end
 re2(wnumob,1)=sqrt(h2((wnumob-1),1))*epsilon((wnumob-1),2);

 re=[re1,re2];
 re=re(numob*0.2+1:wnumob,:);
end