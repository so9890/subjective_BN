%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: data.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes empirical values of moments used for calibration. Saves
% moments in momCheck.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
home = 0; 
if home ==1
    bpath = '/Users/steph/Dropbox/Research/Price_shocks/Matlab_files/Uncertainty/AejFiles';
else
    bpath = '/Users/sfried/Dropbox/Research/Price_shocks/Matlab_files/Uncertainty/AejFiles';
end

cd(bpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = xlsread('matlab4Energy.xls');
pCoal = N(:,2); pGas = N(:,3); pOil = N(:, 4); 
racOilImports = N(:, 5); racOilDomestic = N(:,6); cCoal = N(:,7); cGas = N(:, 8); cOil =N(:,9); 
prodCoal = N(:,10); prodGas = N(:,11); prodOil = N(:,12); prodPetrol = N(:,13); importsPet = N(:,14); 

M = xlsread('matlab4.xls'); %GREEN ENERGY INCLUDES NUCLEAR
year = M(:,1); fossilComp = M(:,2); greenComp = M(:, 3); 
nucComp = M(:,4); energyComp = M(:, 5); totalComp = M(:,6); fossilGov = M(:, 7); greenGov = M(:,8); 
nucGov = M(:,9);  energyGov = M(:, 10); totalGov = M(:, 11);  workforce = M(:, 12); fossilLabor = M(:, 13);
gdp = M(:, 14); XfN = M(:, 15);  XN = M(:, 16); deflate = M(:, 17); RGDPPC = M(:,18); petrol = M(:,19); 


%XfN nominal fixed assets in fossil energy
%XN nominal fixed assets in non-energy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct annual moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


total = totalComp + totalGov;
fossilTotal = fossilGov + fossilComp; 
greenTotal = greenGov + greenComp; 
energyTotal =fossilTotal+ greenTotal; 

nucTotal = nucGov + nucComp;
nonTotal = total - fossilTotal - greenTotal;
pOilReal = pOil./deflate*100; 
gdpReal = gdp./deflate.*100; 
racOilImportsReal = racOilImports./deflate*100; 

growGDP = RGDPPC(2:end)./RGDPPC(1:end-1) -1; 

fracSfComp = fossilComp./totalComp; 
fracSgComp = greenComp./totalComp; 
fracSgTot2 = (greenTotal - nucTotal)./total;
fracXf = XfN./XN;
fracLf = fossilLabor./workforce; 
fracNuc = nucTotal./total; 
fracNucG = nucTotal./greenTotal;
fracSfEng = fossilTotal./energyTotal; 
fracSgEng = greenTotal./energyTotal;
fracPetrol = petrol./total; %fraction of total energy expenditures on petroleum products

%Total energy R&D expenditures in 1972 = 750. Split between fossil and
%green not available. Therefore assume that the fraction is the same as it
%was in 1973. 
energyTotal(24) = 750; 
fracSfEng(24) = fracSfEng(25); 
fracSgEng(24) = fracSgEng(25); 
fossilTotal(24)= energyTotal(24)*fracSfEng(24); 
greenTotal(24) = energyTotal(24)*fracSgEng(24); 
fracSfTot = fossilTotal./total; 
fracSgTot = greenTotal./total;  

oilImports = cOil - prodPetrol; 
shareProd =  (pOil.*prodPetrol + pCoal.*prodCoal + pGas.*prodGas).*1e9./(gdp*1e9);
shareConsume =  (pOil.*cOil + pCoal.*cCoal + pGas.*cGas).*1e9./(gdp*1e9); 
shareImports = importsPet.*1e3./(gdp*1e9); 

pf =(pOil.*prodPetrol + pCoal.*prodCoal + pGas.*prodGas)./(prodPetrol + prodCoal + prodGas);
pfStarReal = [pOilReal(1:19); racOilImportsReal(20:end)];
pfStar = [pOil(1:19);   racOilImports(20:end)];
pfReal = pf./deflate*100; 
theta = pfStar./pf;

%NO NUC
fracSgTotNN = (greenTotal - nucTotal)./(total - nucTotal); 
fracSfTotNN = fossilTotal./(total - nucTotal); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct BGP moments (1961-1970) and shock period momnets (1971-1975)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Non-Science
indBGPS = find(year==1961); 
indBGPE = find(year == 1970); 
indShockS = find(year == 1971); 
indShockE = find(year == 1975); 

%Uncomment lines 112-115 to calculate data for the elasticity of innovation 
%for 2003 oil shock robustness (Table 5).Add 'fracSfTot' and fracSgTot' to
%vars list for this calculation. 
%indBGPS = find(year==1998); 
%indBGPE = find(year == 2002); 
%indShockS = find(year == 2003); 
%indShockE = find(year == 2007); 

vars ={'pfReal', 'fracXf', 'fracLf', 'shareProd', 'pfStar',  'shareImports', 'pfStarReal', 'theta', 'gdpReal', 'RGDPPC','year'};

for i = 1: length(vars)
     temp = zeros(2,1);
     s1= strcat('varNow', '=', vars{i}, ';');
     eval(s1); 
     ind1 = indBGPS; 
     ind2 = indBGPE; %BGP period length
     for j = 1:2
         temp(j) =  mean(varNow(ind1: ind2));
         ind1 = indShockS; 
         ind2 =indShockE; 
     end 
     s2= strcat(vars{i}, '5', '=', 'temp', ';');
     eval(s2);
end 

%Replace mean of RGDPPC on BGP with mean over (1966-1970)
RGDPPC5(1) = mean(RGDPPC(indBGPE -4: indBGPE));


svars ={'fracSfComp', 'fracSgComp', 'fracSfTot', 'fracSgTot', 'fracNuc', 'fracNucG', 'fracSgTotNN', 'fracSfTotNN'};

indStartS = find(year== 1972); 
indEndS  = find(year == 1975); 

for i = 1: length(svars)
     temp = zeros(2,1);
     s1= strcat('varNow', '=', svars{i}, ';');
     eval(s1); 
     for j = 2:2
         temp(j) =  mean(varNow(indStartS: indEndS));
     end 
     s2= strcat(svars{i}, '5', '=', 'temp', ';');
     eval(s2);
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile data used for standard error calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yearBGP = year(indBGPS:indBGPE); 
yearS = year(indShockS:indShockE); 
sgS = fracSgTot(indShockS:indShockE); 
sfS = fracSfTot(indShockS:indShockE); 
fossilBGP = shareProd(indBGPS:indBGPE); 
fossilS = shareProd(indShockS: indShockE); 
importsBGP = shareImports(indBGPS:indBGPE); 
importsS = shareImports(indShockS:indShockE); 
petrolBGP = fracPetrol(indBGPS:indBGPE); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('momCheck','shareProd5', 'shareImports5', 'fracSfTot5', 'fracSgTot5', 'pfStarReal5', 'theta', 'fracXf5', 'fracLf5', 'RGDPPC5');  
save('momCheck','yearBGP', 'yearS', 'sgS', 'sfS', 'fossilBGP', 'fossilS', 'importsBGP', 'importsS', 'petrolBGP', '-append');  


disp('momCheck Saved!');

