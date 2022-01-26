%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour



%% include path to package
clc, clear
folder='/home/sonja/Documents/projects/Overconsumption/codding_model/Own/tools'
addpath(genpath(folder))

cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')

% create folder to save figures
%if CHECK IF DOES OR DOES NOT EXIST 
mkdir('figures/Rep_agent');

%% read in model file and variables
model_rep_agent;

%% read in and check steady state file


