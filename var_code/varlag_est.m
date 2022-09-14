

function[beta,residuals,omega, BIC,...
    AIC, LR, p, q, opt_lag,AR_3d]= pick_lag(data, no_of_lags,k)

%**************************************************************************
% PURPOSE: getting optimal lags and greek letter (k,k) matrices and cov_matrix  
%--------------------------------------------------------------------------
% INPUT:
% - data
% - no_of_lags: maxlags first used to determine opt_lags and then replaced
% late with the optimal lags 
% - k :# of variables 

%--------------------------------------------------------------------------
% OUTPUT:
% - BETA: matrix of coefficients matrix(1+k*lags,k)
% - omega: covariance matrix used in the short run identification 
% -  beta_3d: (k,k,lags) matrices
% COMMENT: remember to accumulate whenever you have data that is not in level form
% otherwise we good and this will be used in your plot 
%--------------------------------------------------------------------------

%global k
%Initialize
AIC = zeros(no_of_lags,1);
BIC = zeros(no_of_lags,1);
HQIC = zeros(no_of_lags,1);
omega_col=zeros(no_of_lags,1);

%%
for lags=1:no_of_lags
    Y_tminus1=lagmatrix(data,1:1:lags);
    Y=data(no_of_lags+1:end,:);
    Y_tminus1=Y_tminus1(no_of_lags+1:end,:);

    [T k]=size(Y);

    constant = ones(T,1);
    X = [constant Y_tminus1];

%% OLS regression 
    beta= X\Y;
    fitted =X*beta;
    residuals = Y-fitted;


    
    omega = cov(residuals);
   
    omega_det = det(omega);
    omega_col(lags,1)=log(abs(det(omega)));

    kappa = k^2*lags+k;

    AIC(lags,1) = log(abs(omega_det))+2*kappa/T;
    BIC(lags,1) = log(abs(omega_det))+kappa*log(T)/T;
   % HQIC(lags,1)= log(abs(omega_det))+2*kappa*log(log(T))/T;
end 

[~,p] = min(BIC);
[~,q]= min(AIC);
%[~,r]=min(HQIC);

%% Likelihood Ratio Test 
m=k^2;
LR=0;
opt_lag=0;


for j=0:no_of_lags-2
    omega_nod=omega_col(no_of_lags-j,1);
    omega_one=omega_col(no_of_lags-(j+1),1);
    LR = abs(T*(omega_nod-omega_one));
    if  (1-chi2cdf(LR,m))<.025
        opt_lag = no_of_lags-j;
        break
    else 
        opt_lag=0;
    end 
end 

b = beta(2:end,:)';
AR_3d=NaN*zeros(k,k,no_of_lags);
for i=1:no_of_lags
    AR_3d(:,:,i)=b(:,k*(i-1)+1:k*i);
end


