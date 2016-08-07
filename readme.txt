% This folder contains the behavioral datasets and Matlab code used for analysis of behavioral data from the paper 
% "Distributed patterns of occipito-parietal functional connectivity predict the precision and variability of visual working memory"
% by Elena Galeano Weber, Tim Hahn, Kirsten Hilger, and Christian J. Fiebach (submitted/under revision; citation to be updated)

% Fiebach Cognitive Neuroscience Lab, Department of Psychology, Goethe University Frankfurt, Germany
% August 2016


% The Matlab code has been previously published by Ronald van den Berg to accompany the
% the paper "Factorial comparison of working memory models" by Van den Berg, Awh, and Ma (Psychological Review, 2014)
% and can be also downloaded from http://www.ronaldvandenberg.org/code.html 
% It is replicated here with the kind permission of the original author, Ronald van den Berg.



%%% file description %%%

see http://www.ronaldvandenberg.org/code.html for detailed file descriptions

fit_factorial_model.m	: written by Ronald van den Berg, Oct 2015, for the tutorial "Modeling delayed-estimation data" given at the Sparks Workshop on Active Perceptual Memory.

code_to_fit_models.m	: code by the authors used to fit the EPA, EPF, VPA, and VPF models to their behavioral datasets; computes AICs for each model, statistics of model fits and AICs.




%%% data structure %%%

% datasets.mat contains the cell array "datasets" with 22 cells (i.e., for 22 subjects).

% Each cell contains three fields (n, condition, and errors) and 477 observations (i.e., trials).

% The fields "n" and "condition" contain the set size (i.e., load 1, 3, or 5) on each trial.
% The field "errors" contains the errors on each trial, i.e., the difference between presented and reported color, 
% a continuous measure of response error in degrees, range [-180, 180].



