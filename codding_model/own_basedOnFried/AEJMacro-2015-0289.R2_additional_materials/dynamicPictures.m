%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: Dynamic pictures
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Produce Figure 1 in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
home =0;

if home ==1
    bpath = '/Users/steph/Dropbox/Research/Price_shocks/Matlab_files';
    f_path = '/Users/steph/Dropbox/Research/Price_shocks/JMP/Paper28';
else
    bpath = '/Users/sfried/Dropbox/Research/Price_shocks/Matlab_files';
    f_path = '/Users/sfried/Dropbox/Research/Price_shocks/JMP/Paper28'; 
end

pic =0;

cd(strcat(bpath, '/Uncertainty/AEJFiles'))
load(strcat(bpath, '/Uncertainty/AEJFiles/','baseResults'));

set(0,'DefaultAxesLineStyleOrder','-|--|-.|:'...
    ,'DefaultLineLineWidth', 4, 'DefaultAxesLineWidth', 1', 'DefaultAxesFontSize', 25, ...
    'DefaultAxesFontName', 'Times');


priceRatio = pg(:,2)./pf(:,2);
laborRatio = lg(:,2)./lf(:,2);
techRatio = ag(:,2)./af(:,2);

priceRatioBase = pg(:,1)./pf(:,1);
laborRatioBase = lg(:,1)./lf(:,1);
techRatioBase = ag(:,1)./af(:,1);
n1 = 19;
year5= zeros(n1,1);
year5(1) = 2010;
for i= 2:1:n1
    year5(i) = year5(i-1)+5;
end

figure(1)
plot([2000; 2005; year5(1:n1)],[priceRatio(1); priceRatio(1); priceRatio(1:n1)], ...
    [2000; 2005; year5(1:n1)],[priceRatioBase(1); priceRatioBase(1); priceRatioBase(1:n1)],'--')
title('Price Effect')
legend('Carbon Tax', 'Baseline', 'Location', 'SW')
xlabel('Year')
ylabel('Pg/Pf')
xlim([year5(1)-10, year5(end)]);
if pic ==1
    saveas(gcf,[f_path ,filesep, 'priceEffect'],'epsc');
end

figure(2)
plot([2000; 2005; year5(1:n1)],[laborRatio(1); laborRatio(1); laborRatio(1:n1)], ...
    [2000; 2005; year5(1:n1)],[laborRatioBase(1); laborRatioBase(1); laborRatioBase(1:n1)],'--')
title('Market Size Effect')
%legend('Baseline', 'Tax', 'Location', 'NW')
xlabel('Year')
ylabel('Lg/Lf')
xlim([year5(1)-10, year5(end)]);
if pic ==1
    saveas(gcf,[f_path ,filesep, 'marketSizeEffect'],'epsc');
end


figure(3)
plot([2000; 2005; year5(1:n1)],[techRatio(1); techRatio(1); techRatio(1:n1)], ...
    [2000; 2005; year5(1:n1)],[techRatioBase(1); techRatioBase(1); techRatioBase(1:n1)],'--')
title('Relative Technology')
%legend('Tax', 'Baseline', 'Location', 'NW')
xlabel('Year')
ylabel('Ag/Af')
xlim([year5(1)-10, year5(end)]);
if pic ==1
    saveas(gcf,[f_path ,filesep, 'tech'],'eps');
end


figure(4)
plot([2000; 2005; year5(1:n1)], [sf(1,1)./S*100; sf(1,1)./S*100; sf(1:n1,2)./S*100]...
    , [2000; 2005; year5(1:n1)], [sf(1,1)./S*100; sf(1,1)./S*100; sf(1:n1,1)./S*100], '--')
title('Percent of Fossil Scientists')
xlabel('Year')
ylabel('Percent')
ylim([0, 2.5]);
xlim([year5(1)-10, year5(end)]);
if pic ==1
    saveas(gcf,[f_path ,filesep, 'sf'],'epsc');
end


figure(5)
plot([2000; 2005; year5(1:n1)], [sg(1,1)./S*100; sg(1,1)./S*100; sg(1:n1,2)./S*100]...
    , [2000; 2005; year5(1:n1)], [sg(1,1)./S*100; sg(1,1)./S*100; sg(1:n1,1)./S*100], '--')
title('Percent of Green Scientists')
xlabel('Year')
ylabel('Percent')
xlim([2000, year5(end)]);
ylim([0, 2.5]);
if pic ==1
    saveas(gcf,[f_path ,filesep, 'sg'],'epsc');
end



