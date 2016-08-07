% function p_error = p_error_factorial(pars, error_vec, N, J2k)
%
% Returns predicted error probabilities under the VP-F
% model, given the set size and parameter values specified in pars.
%
% To get predicitions for EP-A, set tau=0, K=0
% To get predicitions for EP-F, set tau=0
% To get predicitions for VP-A, set K=0
%
% INPUT
%   pars           : model parameters: [J1/J1bar, power, kappa_r, tau, K]
%   error_vec      : vector with response errors
%   N              : set size of the trials from which these response errors came
%
% OUTPUT
%  p_error : probabilities of the errors in error_vec under the specified model
%
% Note that the 2nd model factor is not relevant here,

% Written by RvdB, Oct 2015, for the tutorial "Modeling delayed-estimation
% data" given at the Sparks Workshop on Active Perceptual Memory. Please
% report any bugs or comments to ronald.vandenberg@psyk.uu.se.

function p_error = p_error_factorial(pars, error_vec, N, J2k)

% unpack parameters
J1bar = pars(1);   % mean of gamma distribution over precision
power = pars(2);   
kappa_r = pars(3);
tau = pars(4);     % controls the amount of variability in precision (0 = no variability)
K = pars(5);       % limit on number of items that can be remembered (K>=8 = all items remembered)

% step 0: check if we have any slots at all; if not, return guessing probs
if K==0
    p_error = ones(size(error_vec))*1/2/pi;
    return
end

% step 1: compute #items in memory and p(probed item was in memory)
M = min(N,K);       % number of items in memory: N if there are fewer items than slots; otherwise K
p_mem = min(K/N,1);

% step 2: compute p_error for case that probed item was not in memory (guessing)
p_error_not_remembered = 1/(2*pi);

% step 3: compute p_error for case that probed item was in memory
if tau == 0
    % EQUAL PRECISION
    J = J1bar * M^power; % precision of each item in memory 
    kappa_vec = interp1(J2k.J,J2k.kappa,J);
else
    % VARIABLE PRECISION
    % discretize gamma distribution over J into bins with equal probability masses
    n_gamma_bins  = 50;
    X = linspace(0,1,n_gamma_bins+1);
    X = X(2:end)-diff(X(1:2))/2;
    Jbar = J1bar * M^power;  % mean precision of each item in memory
    J_vec = gaminv(X,Jbar/tau,tau);
    kappa_vec = interp1(J2k.J,J2k.kappa,J_vec,'linear','extrap');
end

% compute probability of target estimation error under every kappa
kappa_vec = min(kappa_vec,1e6);   % HACK - but essential, otherwise numerical problems with EP-X models..
kappa_c = sqrt(bsxfun(@plus,kappa_vec.^2+kappa_r^2,bsxfun(@times,2*kappa_vec*kappa_r,cos(error_vec'))));
p_error_mat = bsxfun(@rdivide,besseli0_fast(kappa_c,1),2*pi*besseli0_fast(kappa_vec,1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus,kappa_c,kappa_vec+kappa_r));
p_error_remembered = mean(p_error_mat,2);

% step 4: combine predictions for "item was in memory" and "item was not in memory"
p_error = p_mem*p_error_remembered + (1-p_mem)*p_error_not_remembered;

