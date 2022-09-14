function[IRF]=irf_est(phi,hmax,omega,shock,k)
%**************************************************************************
% PURPOSE: Computing structural  impulse response of the VAR
%--------------------------------------------------------------------------
% INPUT:
% - K: total number of endogenous variables in the VAR
% - shock: identify variable that to which other variables are responding 
% -  PHI: k-dim matrix (K x K x (N+1)) containing dynamic multipliers of
% the VAR
% - hmax: forecast horizon
% - omega: covariance matrix of error terms 
%--------------------------------------------------------------------------
% OUTPUT:
% - IRF: irf coeffiecient matrix (hmax+1,k) the structural response
% variables for each variable at a certain horizon
% COMMENT: remember to accumulate whenever you have data that is not in level form
% otherwise we good and this will be used in your plot 
%--------------------------------------------------------------------------
b_node=chol(omega,"lower");
IRF=zeros(hmax+1,k)
%%
for h=1:hmax+1

    IRF(h,:) = (phi(:,:,h)*b_node*shock)';
    
end 
  
for h=2:hmax
    IRF(h,1)=IRF(h-1,1)+IRF(h,1);
end 
