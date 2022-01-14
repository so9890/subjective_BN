%*********************************
%PROGRAM PLOTTTAU
%**********************************
clc
clear all;

%% SCATTERPLOT
data     = xlsread('tau_US_plotdata.xlsx','A2:C101');
regconst =  1.71666;
regconst_adj = regconst*1.0955; %adjusted to match g expenditures 0.189*Y
regbeta  =  0.8192254;
lpregovinc    = data(:,2);
ldispinc      = data(:,3);
lpregovinc_vec = linspace(7.2,13.4,10000);
ldispinchat    = regconst + regbeta*lpregovinc_vec;

figure(1);
sizedot = 80;
scatter(lpregovinc,ldispinc,sizedot,'b','filled');
hold on
plot(lpregovinc_vec,ldispinchat,'r','Linewidth',2.5);
grid on;
xlabel('Log of Pre-government Income', 'fontsize', 14);
ylabel('Log of Disposable Income', 'fontsize', 14);
xlim([7.2,13.5]); ylim([7.2,13.5]);


%% PLOT MTR AND ATR IMPLIED BY TAU_US
tauvec    = [1-regbeta];
lambdavec = [exp(regconst_adj)];
ybar = 60000;

for j=1:1;
    tau_temp = tauvec(j);
    lambda_temp =lambdavec(j);
end;
for i=1:2100;
        yv(i)=ybar*0.004*i;
        mtax(i) = 1 - lambda_temp*(1-tau_temp)*yv(i)^(-tau_temp);
        atax(i) = 1 - lambda_temp*yv(i)^(-tau_temp);

	    %yv(i) = yv(i)/ybar;
end;

figure(2);
plot(yv,mtax,'r--',yv,atax,'b-','Linewidth',2.5);
%legend('US Marginal (\tau^{US} = 0.183)', 'US Average (\tau^{US} = 0.183)','Location','SouthEast', 'fontsize', 16);
legend('Marginal Tax Rate', 'Average Tax Rate','Location','SouthEast', 'fontsize', 16);
grid on;
ylim([-0.31 0.605]);xlim([0.0 500000]);
xlabel('Pre-government Income', 'fontsize', 14);
ylabel('Tax Rates', 'fontsize', 14);



