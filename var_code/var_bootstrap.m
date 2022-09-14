function[gdp,inf, eps ]=bts(data,lags,k,shock,hmax)
%**************************************************************************
% PURPOSE: Bootstrapping to create confidence intervals 
%--------------------------------------------------------------------------
% INPUT:
% - data: variables that we are interested in
% - shock: identify variable that to which other variables are responding 
% -  lags: lags, or q 
% - hmax: forecast horizon
% - k: # of variables  
%--------------------------------------------------------------------------
% OUTPUT: var1, var2, ...,vark
% depends on the number of variables 
% variables (hmax+1,N_rep) matrix 
%--------------------------------------------------------------------------

r=length(data);
[beta1,res]= lag_est(data,lags,k);
T_burnin=100;
%%
eps=[];
gdp=[];
inf=[];

c_hat=beta1(1,:);
u_hat=res(lags+1:end,:);
psi =beta1(2:end,:);

%%
for m=1:500
    u_sample =datasample(u_hat,r+T_burnin);
    possible_values= lags:length(u_hat);
    pos = randi(length(possible_values));
    tau = possible_values(1,pos);
    init_cond = flipud(data(tau-lags+1:tau,:));

    y_init=init_cond(1,:);

%%
    for i=2:lags
        y_init=[y_init,init_cond(i,:)];
    end
%% 

    new_y=[];
    for j=1:r+T_burnin
        y_gen=y_init*psi +c_hat+u_sample(j,:);
        new_y=[new_y; y_gen];
        y_init = [y_gen y_init(1:end-k)]; %push the observations
    end 
    new_y=new_y(T_burnin+1:end,:);
    
    [boot,res_boot,omega_boot,~,~,~,~,~,~,boot_3d]=lag_est(new_y,lags,k);
    PHI_boot = dyn_multipliers(k,lags,boot_3d,hmax);
    big_theta=irf_est(PHI_boot,hmax,omega_boot,shock,k);
    
    gdp=[gdp big_theta(:,1)];
    inf=[inf big_theta(:,2)];
    eps=[eps big_theta(:,3)];

end    
    
    