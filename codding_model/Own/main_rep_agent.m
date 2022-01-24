%% Representative Agent Model

% with two skill types, 
% higher disutility from labour for high skill labour



%% include path to package
folder='/home/sonja/Documents/projects/Overconsumption/codding_model/Own/tools'
addpath(genpath(folder))

cd('/home/sonja/Documents/projects/Overconsumption/codding_model/Own')

% create folder to save figures
%if CHECK IF DOES OR DOES NOT EXIST 
mkdir('figures/Rep_agent');

%% read in model file and variables
model_rep_agent;

%% read in and check steady state file


% check if ss values set model equations numerically close to zero
cell.y=arrayfun(@char, list.y, 'uniform', 0);
cell.x=arrayfun(@char, list.x, 'uniform', 0);
cell.yp=arrayfun(@char, list.yp, 'uniform', 0);
cell.xp=arrayfun(@char, list.xp, 'uniform', 0);

modelSS=subs(model.f, [cell.x, cell.xp, cell.yp, cell.y], [Xss', Xss', Yss', Yss']);
