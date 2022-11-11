%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: standardErrors.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Calculates bootstrapped standard errors for parameters
% calibrated with the method of moments procedure
% Calls program: calibrate.m
% saves results: se4.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
home = 1;
if home ==1
    bpath = '/Users/steph/Dropbox/Research/Price_shocks/Matlab_files/Uncertainty/AEJFiles';
else
    bpath = '/Users/sfried/Dropbox/Research/Price_shocks/Matlab_files/Uncertainty/AEJFiles';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load momCheck
load base

p0 = [alphag, eta, deltaye, epsf, deltaff, nu1];
deltaPfStarD = (pfStarReal5 - pfStarReal5(1))./pfStarReal5(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate standard deviation of BGP moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdF = std(fossilBGP);  %standard deviation of fossil share along BGP
sdI = std(importsBGP);
sdP = std(petrolBGP); %standard deviation of petroleum R&D as a share of the total on the BGP

sdPF = sdP*fracSfTot5(2)/mean(petrolBGP); %Rescale sdP to account for different means between fossil share of R&D and petroleum share of R&D
sdPG = sdP*fracSgTot5(2)/mean(petrolBGP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resample moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reps = 300;
momF1 = zeros(reps,1); momI1 = zeros(reps,1);
momF2 = zeros(reps,1); momI2 = zeros(reps,1);
momSf = zeros(reps,1); momSg = zeros(reps,1);

for i = 1:reps
    momF1(i) = shareProd5(1) + normrnd(0, sdF); %the observed BGP mean plus the standard deviation of the data in the BGP
    momI1(i) = shareImports5(1) + normrnd(0, sdI); %the observed BGP mean plus the standard deviation of the data in the BGP
    momF2(i) = shareProd5(2) + normrnd(0, sdF); %the observed shock period mean plus the standard deviation of the data in the BGP
    momI2(i) = shareImports5(2) + normrnd(0, sdI); %the observed shock period mean plus the standard deviation of the data in the BGP
    momSf(i) = fracSfTot5(2)  + normrnd(0, sdPF); %the observed shock period mean plus the standard devaition of the petrol data in the BGP
    momSg(i) = fracSgTot5(2)  + normrnd(0, sdPG); %the observed shock period mean plus the standard devaition of the petrol data in the BGP
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recalibrate parameters for each resampled moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = zeros(reps, 7);
dist = zeros(reps, 1);

for i = 1:reps
    momD =[momI1(i); momI2(i); momF1(i); momF2(i); momSf(i); momSg(i)]*100;
    
    [dist(i), params(i, :)] = calibrate(momD, p0,alphaf, alpham, deltaef, epse, epsf, epsy, L,phi ,rhonf, rhong,...
        S, theta0, deltaPfStarD, guess0, guessBgp, options, Ttarg, i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate standard errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Throw out runs that did not converge
mat = [dist, params];
matSort = sortrows(mat, 1);
tooBig = find(matSort(:,1) > 1e-15, 1, 'first');
matSort2 = matSort(1:tooBig -1, :);
%deltayeT = (matSort2(:,4)./(1-matSort2(:,4))).^epsy; %transformed deltaye
se = std(matSort(:, 2:end));

%[alphag, eta, deltaye, epsf, deltaff, nu1, gamma];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('standardErrors', 'alphag', 'eta', 'deltaye', 'epsf', 'deltaff', 'nu1', 'gamma', 'se');
save se4;



