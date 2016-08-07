% Code for fitting the EPA, EPF, VPA, and VPF model to the behavioral datasets of the paper 
% "Distributed patterns of occipito-parietal functional connectivity predict the precision of visual working memory" 
% by Elena Galeano Weber, Tim Hahn, Kirsten Hilger, and Christian J. Fiebach (submitted/under revision; citation to be updated)

% Fiebach Cognitive Neuroscience Lab, Department of Psychology, Goethe University Frankfurt, Germany
% August 2016


% The Matlab functions embedded in this code have been previously published by Ronald van den Berg to accompany 
% the paper "Factorial comparison of working memory models" by Van den Berg, Awh, and Ma (Psychological Review, 2014)
% and can be also downloaded from http://www.ronaldvandenberg.org/code.html 
% They are replicated here with the kind permission of the original author, Ronald van den Berg.


clear all
close all

%% load data
load('datasets.mat');

%% transform errors in a range from -pi to pi and define data structure
for ii=1:length(datasets)
    
    data{ii}.error = datasets{1,ii}.errors*pi/180;
    data{ii}.N = datasets{ii}.n
    data{ii}.dist_error_vec = nan;
end

%% fit model
for ii=1:length(datasets)
    
    
    fprintf('Subject #%d\n',ii);
    [fitpars_EPA(ii,:) CI_lower_EPA(ii,:) CI_upper_EPA(ii,:) log_lh_EPA(ii)] = fit_factorial_model([1 1],data{ii});
    [fitpars_EPF(ii,:) CI_lower_EPF(ii,:) CI_upper_EPF(ii,:) log_lh_EPF(ii)] = fit_factorial_model([1 2],data{ii});
    [fitpars_VPA(ii,:) CI_lower_VPA(ii,:) CI_upper_VPA(ii,:) log_lh_VPA(ii)] = fit_factorial_model([2 1],data{ii});
    [fitpars_VPF(ii,:) CI_lower_VPF(ii,:) CI_upper_VPF(ii,:) log_lh_VPF(ii)] = fit_factorial_model([2 2],data{ii});
    
    
end

%% compute AIC values
AIC_EPA = -2*log_lh_EPA + 2*3;
AIC_EPF = -2*log_lh_EPF + 2*4;
AIC_VPA = -2*log_lh_VPA + 2*4;
AIC_VPF = -2*log_lh_VPF + 2*5;

clear ii data 

save('fit_results.mat');



%-----------------------------------------------------------------------------------------------------%
%                                       Descriptive Statistics                                        %
%-----------------------------------------------------------------------------------------------------%
%% compute medians of AIC values
AIC_EPA_med = median(AIC_EPA);
AIC_EPF_med = median(AIC_EPF);
AIC_VPA_med = median(AIC_VPA);
AIC_VPF_med = median(AIC_VPF);

AIC_all_medians = [AIC_EPA_med,AIC_EPF_med,AIC_VPA_med,AIC_VPF_med];

%% Print Table of individual AIC-values
Table_AICs = [AIC_EPA',AIC_EPF',AIC_VPA',AIC_VPF']
Table_AICs_header = {'AIC_EPA','AIC_EPF','AIC_VPA','AIC_VPF'}

%-----------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------%

%% model parameters of 'variable precision fixed capacity' (VPF) model
J1bar = fitpars_VPF(:,1);  % mn precision
power = fitpars_VPF(:,2); % decline in mean precision with load
kappa_r = fitpars_VPF(:,3); % response noise
tau = fitpars_VPF(:,4); % var precision
K = fitpars_VPF(:,5);  % maximum number 


%% relationship between mn precision and set size (ss) is defined in a power-law fashion
for vpInd = 1:length(J1bar)
J1bar_ss1(vpInd,1) = J1bar(vpInd,1)*1.^power(vpInd,1); 
J1bar_ss3(vpInd,1) = J1bar(vpInd,1)*3.^power(vpInd,1); 
J1bar_ss5(vpInd,1) = J1bar(vpInd,1)*5.^power(vpInd,1); 
end

precision=[J1bar_ss1,J1bar_ss3,J1bar_ss5];


%% plot load-dependent decline of mn precision
precision=log(precision);
y = mean(precision);
e = std(precision)/sqrt(length(precision(:,1))); % standard error of the mean
x=[1 3 5]

figure;hold on
errorbar(y,e,'-','Color', 'k',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',6);
                     
set(gca,'XTick',1:1:3)
set(gca,'XTickLabel',{'load 1','load 3','load 5'}, 'FontSize', 16)
ylabel('log(mean precision)','FontSize',16)


%% logtransform mn precision, var precision, and response noise
J1bar = log(J1bar);
tau = log(tau);
kappa_r = log(kappa_r);

%% Print Table of individual model parameters of the VPF model [K, (log(mn precision, power, 
%% log(var precision), log(response noise)]
VPFmodel_params_all = [K, J1bar, power, tau, kappa_r]
VPFmodel_params_all_header = {'K', 'log(mn precision)', 'power', 'log(var precision)', 'response noise'}

%-----------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------%

%% Calc median, mean, sem, sd, min., max.
% mn precision
J1bar_mn = mean(J1bar);
J1bar_med = median(J1bar);
J1bar_sem = std(J1bar)/sqrt(length(J1bar)); % standard error of the mean
J1bar_sd = std(J1bar);
J1bar_min = min(J1bar);
J1bar_max = max(J1bar);
J1bar_des = [J1bar_med,J1bar_mn,J1bar_sem,J1bar_sd,J1bar_min,J1bar_max];

% power
power_mn = mean(power);
power_med = median(power);
power_sem = std(power)/sqrt(length(power)); 
power_sd = std(power);
power_min = min(power);
power_max = max(power);
power_des = [power_med,power_mn,power_sem,power_sd,power_min,power_max];

% kappa_r/response noise
kappa_r_mn = mean(kappa_r);
kappa_r_med = median(kappa_r);
kappa_r_sem = std(kappa_r)/sqrt(length(kappa_r)); 
kappa_r_sd = std(kappa_r);
kappa_r_min = min(kappa_r);
kappa_r_max = max(kappa_r);
kappa_r_des = [kappa_r_med,kappa_r_mn,kappa_r_sem,kappa_r_sd,kappa_r_min,kappa_r_max];

% tau/var precision
tau_mn = mean(tau);
tau_med = median(tau);
tau_sem = std(tau)/sqrt(length(tau)); 
tau_sd = std(tau);
tau_min = min(tau);
tau_max = max(tau);
tau_des = [tau_med,tau_mn,tau_sem,tau_sd,tau_min,tau_max];

% K
K_mn = mean(K);
K_med = median(K);
K_sem = std(K)/sqrt(length(K)); 
K_sd = std(K);
K_min = min(K);
K_max = max(K);
K_des = [K_med,K_mn,K_sem,K_sd,K_min,K_max];


%% Print Table 1 %%
Table_01_pars = [K_des;J1bar_des;power_des;tau_des;kappa_r_des]
Table_01_header_col = {'median', 'mean', 'SEM', 'SD', 'Min.', 'Max.'} % header of columns
Table_01_header_row = {'K';'log(mn precision)'; 'power'; 'log(var precision)';'log(response noise)'} % header of rows


%-----------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------%


% %% plot model fit
% vpInd = 1; % number of subject (ranges from 1 to 22)
% fitpars = fitpars_VPF(vpInd,:);
% plot_error_histograms(data{1,vpInd},fitpars);






