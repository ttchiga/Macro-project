%%% VAR SCRIPT
%ii)
PI= diff(PCE);
VAR_IP= IP(2:end,1);
MPv_upd=MP_upd(2:end,1);
Y_t=[VAR_IP PI MPv_upd]; 
[T k]=size(Y_t);
shock=zeros(k,1);
shock(3,1)=1;
[beta,residuals,~, ~,~, ~, p2, q2]= varlag_est(Y_t, 12,k);
%take q2
%%
%iii)


[beta,residuals,omega, ~,~, ~, ~, ~, ~,beta_3d]= varlag_est(Y_t, q2,k);
PHI = var_multipliers(k,q2,beta_3d,hmax);
IRF_coeff=varirf_est(PHI,hmax,omega,shock,k);

%%
[gdp, inflation, mp ]=var_bootstrap(Y_t,q2,k,shock,hmax);
errband_gdp = prctile(gdp,[5  95],2);
errband_inf = prctile(inflation,[5  95],2);
errband_eps = prctile(mp,[5  95],2);

%% PLOT Bootstrap CI
%IRFandCIplot(big_theta(:,1),errband_tax(:,1),errband_tax(:,2),'IRF_bootstrap_IP.png',12);
IRFandCIplot(IRF_coeff(:,1),errband_gdp(:,1),errband_gdp(:,2),'IRF_bootstrap_gdp.png',hmax);
IRFandCIplot(IRF_coeff(:,2),errband_inf(:,1),errband_inf(:,2),'IRF_bootstrap_inf.png',hmax);
IRFandCIplot(IRF_coeff(:,3),errband_eps(:,1),errband_eps(:,2),'IRF_bootstrap_mp.png',hmax);