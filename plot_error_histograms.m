% function plot_error_histograms(data, parameters)
% 
% Plots the error histograms of a given data set, separately for each set size.
% If parameters are provided, then a model prediction is plotted on top of
% the histogram.
%
% INPUT
%  data.N        : vector with set size for each trial
%  data.error    : vector with response error for each trial
%  parameters    : vector with parameter values, [kappa1bar, tau, K] (optional)
%
% The parameter contains all parameters of the most general model, VP-F.
% To plot EP-A predictions, set tau=0 and K=8
% To plot EP-F predictions, set tau=0
% To plot VP-A predictions, set K=8

% Written by RvdB, Oct 2015, for the tutorial "Modeling delayed-estimation 
% data" given at the Sparks Workshop on Active Perceptual Memory. Please 
% report any bugs or comments to ronald.vandenberg@psyk.uu.se.

function plot_error_histograms(data, parameters)

if exist('parameters','var')
    kappa1bar = parameters(1);
    tau = parameters(2);
    K = parameters(3);
end

% open a large figure window
figure
set(gcf,'Position',get(gcf,'Position').*[0.3 0.3 1.5 0.75],'PaperPosition',get(gcf,'PaperPosition').*[0.3 0.3 1.5 0.75]);

% loop over set sizes to plot histograms and compute circ var and circ kurt
maxY = -Inf;
uN = unique(data.N);
J2k.kappa = [linspace(0,10,250) linspace(10.001,1e4,250)];
J2k.J = J2k.kappa.*besseli(1,J2k.kappa,1)./besseli(0,J2k.kappa,1);
for ii=1:numel(uN)
    % plot histogram
    subplot(1,3,ii);
    set(gca,'FontSize',16);
    hold on
    X = linspace(-pi,pi,25);
    X = X(2:end)-diff(X(1:2))/2;
    Y_emp = hist(data.error(data.N==uN(ii)), X);
    Y_emp = Y_emp/sum(Y_emp)/diff(X(1:2));
    bar(X,Y_emp,'k');
    X = linspace(-pi,pi,250);
    X = X(2:end)-diff(X(1:2))/2;
    Y_fit = p_error_factorial(parameters,X,uN(ii),J2k);
    plot(X,Y_fit,'r','linewidth',3);
    xlim([-pi pi]);
%     maxY = max(maxY,max([Y_fit, Y_emp]));
    maxY = max(maxY,max(Y_emp));
end

for ii=1:3
    subplot(1,3,ii);
    ylim([0 maxY*1.3]);
    xlim([-pi pi]);
    set(gca,'Xtick',[-pi 0 pi],'XtickLabel',{'-pi','0','pi'});
%     title(sprintf('Set size %d',ii));
    %text(-3,maxY*1.2,sprintf('Set size %d',ii),'fontsize',10);
    if ii>4
        xlabel('response error');
    else
        set(gca,'Xtick',[]);
    end
    if ii==1 || ii==5
        ylabel('probability','fontsize',16);
    else
        set(gca,'Ytick',[],'fontsize',16);
    end
end