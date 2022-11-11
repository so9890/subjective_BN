%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: cev.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: computes the difference in welfare between the baseline and the
% cabon tax when consumption in the baseline is scaled by factor lambda. 
% Inputs: 
% Guess for CEV: lambda
% Consumption in baseline: cBase
% Consumption under tax: cTax
% Disount rate (detrended and adjust for five year time period: discount
% Utility function parameter: theta
% Outputs: 
% Welfare difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function answ = cev(lambda, cBase, cTax, discount, theta)

%Start CEV calculation in 2015

utils = (lambda*cBase).^(1-theta)./(1-theta); 
utilsTax = (cTax).^(1-theta)./(1-theta); 

welfareTax =0;
welfare =0; 
for i = 1:length(cTax)
    welfareTax = welfareTax + discount^i*utilsTax(i); 
    welfare = welfare + discount^i*utils(i);
end 

answ =welfare -welfareTax;
