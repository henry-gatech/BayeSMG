%%% This function executes the BayeSMG matrix completion
function [X_hat,lb,ub] = BayeSMG(X_omega,omega,r,eta)
%%%
% X_omega: m1*m2 matrix containing noisy observations of matrix entries
% omega: m1*m2 matrix containing the indices of observed entries, 1/0
% r: scalar, Rank of the matrix
% eta: scalar, observation noise deviation
%%%
% load the sparco package
cd sparco-1.2
run sparcoSetup.m
cd '../'

% setup the observed indices
m1 = size(X_omega,1);
m2 = size(X_omega,2);
[idxx,idxy] = find(omega);
n = sum(omega(:));
obs = zeros(n,1);
for i = 1:n
    obs(i) = X_omega(idxx(i),idxy(i));
end

% setup the hyper-parameters (values can be changed)
eta_flg = true; %is eta known?
as2 = 0.01; %\alpha for \sigma^2 prior (weakly informative)
bs2 = 0.01; %(rate) \beta for \sigma^2 prior (weakly informative)
nsamp = 1000; %number of MCMC samples
nu = 1; %t-dist d.f. for independence sampler
burn_gql = 10; %burnin for independence sampler
burn = nsamp/5; %burnin for bsmg
gibbs_flg = true; %gibbs for manifold sampling
e2i = eta^2; %initial \eta^2
s2i = 1; %initial \sigma^2
hpd_perc = 0.95; %HPD percentage

% Run MCMC
samp_bsmg = mcmc_smg(omega,obs,r,nsamp,true,nu,burn_gql,...
    e2i,s2i,as2,bs2,eta_flg,0,0,gibbs_flg,[]);
% estimate the matrix
pm_bsmg = reshape(mean(samp_bsmg.XX(burn,:,:),1),[m1 m2]);
% acquire the UQ
[lw_smg, up_smg] = hpd(samp_bsmg.XX,burn,hpd_perc);
% setup the outputs
X_hat = pm_bsmg;
lb = lw_smg;
ub = up_smg;

end