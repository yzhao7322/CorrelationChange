%% This package provides the coding programme to implement the tests in the following paper.
% Reference: M. Barassi, L. Horváth & Y. Zhao (2020) Change‐Point Detection in the Conditional Correlation Structure of Multivariate Volatility Models, Journal of Business & Economic Statistics.
% Y.ZHAO, Apr/2018

% The empcv_MGARCH.m function calculates the empirical critical values for M_T(1) and M_T(2) statistics, whose limits are discussed in Theorem 1-2.
% input: nobs - sample size
%        reps - number of replications to obtain the critical values (usually set 2000-5000)
%        model_type - the specification of the underlying parametric model, switch it from
%                     'bekk' - BEKK-GARCH model
%                     'factor' - Factor-GARCH model
%                     'ccc' - CCC-GARCH model
%                     'dcc' - DCC-GARCH model
%                     'cdcc' - CDCC-GARCH model
% output: cv_M1 - critical values of M_T(1) test at 90%, 95%, 99% significance levels
%         cv_M2 - critical values of M_T(2) test at 90%, 95%, 99% significance levels
% note: the longrun covariance is computed with the Bartlett kernel function, which could be altered.
function [cv_M1, cv_M2]=empcv_MGARCH(nobs,reps,model_type)
cv_M1=nan(3,1);cv_M2=nan(3,1);
lctr=nan(reps,1);octr=nan(reps,1);

switch lower(model_type)
    case{'bekk'}
        parfor i=1:reps
          [re]=concorrsim(0,nobs);
          [lambdab, omegab, ~, ~] = NorCUSUM(re,'bekk','bartlett');
          lctr(i,1)=lambdab;
          octr(i,1)=omegab;
        end
    case{'factor'}
        [parfor i=1:reps
          [re]=concorrsim(0,nobs);
          [lambdaf, omegaf, ~, ~] = NorCUSUM(re,'factor','bartlett');
          lctr(i,1)=lambdaf;
          octr(i,1)=omegaf;
        end
    case{'ccc'}
        parfor i=1:reps
          [re]=concorrsim(0,nobs);
          [lambdac, omegac, ~, ~] = NorCUSUM(re,'ccc','bartlett');
          lctr(i,1)=lambdac;
          octr(i,1)=omegac;
        end
    case{'dcc'}
        parfor i=1:reps
          [re]=concorrsim(0,nobs);
          [lambdad, omegad, ~, ~] = NorCUSUM(re,'dcc','bartlett');
          lctr(i,1)=lambdad;
          octr(i,1)=omegad;
        end
    case{'cdcc'}
        parfor i=1:reps
          [re]=concorrsim(0,nobs);
          [lambdaa, omegaa, ~, ~] = NorCUSUM(re,'cdcc','bartlett');
          lctr(i,1)=lambdaa;
          octr(i,1)=omegaa;
        end
    otherwise
        error('The model type must be "bekk", or "factor", or "ccc", or "dcc", or "cdcc".');
end

cv_M1(1,1)=quantile(lctr,0.90);cv_M1(2,1)=quantile(lctr,0.95);cv_M1(3,1)=quantile(lctr,0.99);
cv_M2(1,1)=quantile(octr,0.90);cv_M2(2,1)=quantile(octr,0.95);cv_M2(3,1)=quantile(octr,0.99);

end



%% The NorCUSUM.m compute the statistics M_T(1) and M_T(2).
% input: data - the objective data of interest
%        model_type - the specification of the underlying parametric model, switch it from
%                     'bekk' - BEKK-GARCH model
%                     'factor' - Factor-GARCH model
%                     'ccc' - CCC-GARCH model
%                     'dcc' - DCC-GARCH model
%                     'cdcc' - CDCC-GARCH model
%        kernel_type - the kernel function that is chosen to compute the longrun covariance
%                      'bartlett'- Bartlett kernel function
%                      'parzen' - Parzen kernel function
%                      'truncated' - Truncated kernel function
%                      'tukey'- Tukey-Hanning kernel function
%                      'qs' - Quadratic Spectral kernel function 
%                      'flat' - Flat-top kernel function
% output: M1 - M_T(1) statistics
%         M2 - M_T(2) statistics
%         ind - detected change point under the alternative hypothesis
%         mat - the sequence of CUSUM statistics.
function [M1, M2, ind, mat] = NorCUSUM(data,model_type,kernel_type)
[nobs,d]=size(data);
devoldata=zeros(nobs,d);
dm=mean(data);
for i=1:d
    data(:,i)=data(:,i)-dm(1,i);
end
switch lower(model_type)
        // case{'vec'}
        //     [~, ~, ht, ~, ~, ~, ~] = scalar_vt_vech(data,[],1,0,1,'Diagonal');
        case{'bekk'}
            [~, ~, ht, ~, ~] = bekk(data,[],1,0,1,'Full');
        // case{'abekk'}
        //     [~, ~, ht, ~, ~] = bekk(data,[],1,1,1,'Full');
        case{'factor'}
            [~,ht,~,~,~]=o_mvgarch(data,2,1,0,1);
        case{'ccc'}
            [~, ~, ht, ~, ~] = ccc_mvgarch(data,[],1,0,1);
        case{'dcc'}
            [~, ~, ht, ~, ~, ~, ~, ~]=dcc(data,[],1,0,1,1,0,1);
        // case{'adcc'}
        //     [~, ~, ht, ~, ~, ~, ~, ~]=dcc(data,[],1,1,1,1,1,1);
        case{'cdcc'}
            [~, ~, ht, ~, ~, ~, ~, ~]=cdcc(data,[],1,0,1,1,1,1,[],'2-stage');   
        otherwise
            error('Incorrect model_type choice.')
end
for i=1:d
    v=zeros(nobs,1);
    v(:,1)=ht(i,i,:);
    devoldata(:,i)=data(:,i)./sqrt(v);
end
switch lower(kernel_type)
    case{'bartlett'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'bartlett');
    case{'parzen'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'parzen');
    case{'truncated'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'truncated');
    case{'tukey'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'tukey');
    case{'qs'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'qs');
    case{'flat'}
        [M1, M2, ind, mat]=LambdaOmega(devoldata,'flat');
    otherwise
        error('The kernel type must be "Bartlett", or "Parzen", or "Tukey", or "qs", or "tf", or ''flat''.');
end
end

