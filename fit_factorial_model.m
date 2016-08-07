% function [fitpars, CI_low, CI_up, log_lh] = fit_factorial_model(model_flags,data)
%
% Fits a factorial model to delayed-estimation data. 
%
% INPUT
%  model_flags(1) : 1=remembered items have equal precision, 2=variable precision
%  model_flags(2) : 1=all items remembered, 2=fixed slot limit
%  data.N         : vector with set size for each trial
%  data.error     : vector with response error for each trial
% 
% OUTPUT
%  fitpars        : Maximum-likelihood (ML) estimates of the model 
%                   parameters, [J1bar, power, kappa_r, tau, K]. In EP models, returned tau
%                   equals 0 and in -A models, K equal Inf.
%                   
%  log_lh         : log likelihood of the the ML parameter estimates
%

% Written by RvdB, Oct 2015, for the tutorial "Modeling delayed-estimation data" 
% given at the Sparks Workshop on Active Perceptual Memory. Please 
% report any bugs or comments to ronald.vandenberg@psyk.uu.se.

function [fitpars, fitpars_CI_lower, fitpars_CI_upper, log_lh] = fit_factorial_model(model_flags,data)

% Unfortunately, we cannot use fminsearch to find the maximum-likelihood 
% estimate of K, because fminsearch cannot handle integer variables. 
% However, since our maximum set size is 8, we need to consider only 9 
% possible values for K: 0, 1, ..., 8. (a model with K>8 will give the same 
% predictions as K=8, as in all those cases, all items are remembered).
% 
% Therefore, to find the maximum-likelihood estimate of K, we can do the 
% following:
% 1. Loop over all possible values of K, i.e., K=0, 1, ..., 8
% 2. For each K, optimize the remaining parmater(s) using fminsearch
% 3. Determine which value of K maximized the log likelihood and return 
%    those results as the global maximum

factor1_names = {'EP','VP'};
factor2_names = {'A','F'};

% we first build a vector of K-values to consider 
if model_flags(2) == 1;
    K_vec = max(data.N); % all items encoded, which is the same as having a slot limit equal to the maximum set size
elseif model_flags(2) == 2
    K_vec = 0:max(data.N); % we consider all possible values of K
end

% get mapping between kappa and J
J2k.kappa = [linspace(0,10,250) linspace(10.001,1e4,250)];
J2k.J = J2k.kappa.*besseli(1,J2k.kappa,1)./besseli(0,J2k.kappa,1);

% loop over all possible values of K
fprintf('Fitting model %s-%s\n',factor1_names{model_flags(1)},factor2_names{model_flags(2)});
for ii=1:numel(K_vec)    
    % find ML estimates of the remaining parameters, given our current value of K
    if model_flags(1) == 1
        % EP: run fminsearch with initial estimate J1=50, power=-1, kappa_r=100
        [fp, neg_log_lh] = fminsearch(@(pars) -log_lh_function(pars, model_flags, K_vec(ii), data, J2k), [50 -1 100], optimset('Display', 'off'));
        pars_tmp(ii,:) = fp;
        fitpars_all(ii,:) = [fp 0 K_vec(ii)];  % [J1, power, kappa_r, tau, K]
    elseif model_flags(1) == 2
        % VP: run fminsearch with initial estimates J1bar=50, power=-1, kappa_r=100, tau=19
        [fp, neg_log_lh] = fminsearch(@(pars) -log_lh_function(pars, model_flags, K_vec(ii), data, J2k), [50 -1 100 10], optimset('Display', 'off'));
        pars_tmp(ii,:) = fp;
        fitpars_all(ii,:) = [fp, K_vec(ii)]; % [J1bar, power, kappa_r, tau, K]
    end
    % store the log likelihood and parameter estimates 
    log_lh_all(ii) = -neg_log_lh;   
end

% now determine which value of K gave the highest likelihood
best_idx = find(log_lh_all==max(log_lh_all));
fitpars = fitpars_all(best_idx,:);
log_lh = log_lh_all(best_idx);

% get Hessian and compute parameter uncertainty
warning off
[A B C D E Hessian] = fminunc(@(pars) -log_lh_function(pars, model_flags, K_vec(best_idx), data, J2k), pars_tmp(best_idx,:),  optimset('Display', 'off'));
warning on
if any(diag(Hessian)<0)
    fitpars_std = NaN(size(fitpars));
else
    fitpars_std = sqrt(diag(Hessian.^-1))';
end
fitpars_CI_lower = NaN(size(fitpars));
fitpars_CI_upper = NaN(size(fitpars));
fitpars_CI_lower(1:numel(fitpars_std)) = fitpars(1:numel(fitpars_std))-1.96*fitpars_std;
fitpars_CI_upper(1:numel(fitpars_std)) = fitpars(1:numel(fitpars_std))+1.96*fitpars_std;

% return K=Inf in case we fitted EP-A or VP-A
if model_flags(2)==1
    fitpars(end)=Inf;
end

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% This function returns the log likelihood of a set of parameters for a given data set  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
function log_lh = log_lh_function(pars, model_flags, K, data, J2k)

% build parameter vector for p_error_factorial: [J1bar, tau, K]
par_vec(1) = pars(1); % J1 / J1bar
par_vec(2) = pars(2); % power
par_vec(3) = pars(3); % kappa_r
if model_flags(1)==1
    % we are fitting an EP model -> set tau=0
    par_vec(4) = 0;
elseif model_flags(1)==2
    % we are fitting a VP model -> tau is a free parameter, specified as input to this function
    par_vec(4) = pars(4);
end
par_vec(5) = K; % add K to the parameter vector

% return log_lh = -Inf if any of the parameters is negative (all model parameters can only take positive values)
if any(par_vec([1 3 4 5])<0)
    log_lh=-Inf;
    return
end

% compute log likelihood for each trial by looping over set sizes
uN = unique(data.N); % list of unique set sizes in the data set
log_lh = 0; % initialize function output to 0
for ii=1:numel(uN)    
    % set current N
    N = uN(ii);
    
    % get response errors from all trials with set size equal to uN(ii)
    error_vec = data.error(data.N==N);
    
    % compute probability of each error 
    p_error = p_error_factorial(par_vec, error_vec, N, J2k);
      
    % set all p=0 to p=0.00001 (to avoid log(0) below, which would give problems)
    p_error = max(p_error, 1e-5); 
    
    % compute summed log likelihood over all trials with this set size and update total
    log_lh = log_lh + sum(log(p_error));     
end
% fprintf('pars=[%2.1f, %2.2f, %2.1f, %2.1f, %2.1f], log lg=%2.1f\n',par_vec,log_lh);