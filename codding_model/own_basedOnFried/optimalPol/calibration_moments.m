function outputFossil = calibration_moments()
% calibration

% this file 
% 1) reads in data and generates moments 
% 2) matches model equations in baseyear

% CURRENTLY USES FRIED DATA
clear
addpath('data')

%- read in data
dataF = xlsread('matlab4Energy.xls');
pCoal = dataF(:,2); pGas = dataF(:,3); pOil = dataF(:, 4); 
racOilImports = dataF(:, 5); racOilDomestic = dataF(:,6); cCoal = dataF(:,7); cGas = dataF(:, 8); cOil =dataF(:,9); 
prodCoal = dataF(:,10); prodGas = dataF(:,11); prodOil = dataF(:,12); prodPetrol = dataF(:,13); importsPet = dataF(:,14); 

dataF2 = xlsread('matlab4.xls'); %GREEN ENERGY INCLUDES NUCLEAR
year = dataF2(:,1); fossilComp = dataF2(:,2); greenComp = dataF2(:, 3); 
nucComp = dataF2(:,4); energyComp = dataF2(:, 5); totalComp = dataF2(:,6); fossilGov = dataF2(:, 7); greenGov = dataF2(:,8); 
nucGov = dataF2(:,9);  energyGov = dataF2(:, 10); totalGov = dataF2(:, 11);  workforce = dataF2(:, 12); fossilLabor = dataF2(:, 13);
gdp = dataF2(:, 14); XfN = dataF2(:, 15);  XN = dataF2(:, 16); deflate = dataF2(:, 17); RGDPPC = dataF2(:,18); petrol = dataF2(:,19); 

%- calculate moments

%-- share of fossil production/consumption in total GDP
outputFossil= (prodPetrol+prodCoal + prodGas); % in model taken as the same => no weigthing 
shareProd =  (pOil.*prodPetrol + pCoal.*prodCoal + pGas.*prodGas).*1e9./(gdp*1e9);
shareConsume =  (pOil.*cOil + pCoal.*cCoal + pGas.*cGas).*1e9./(gdp*1e9); 

%-- from data on RD estimate green scientists => sg/sf => Ag/Af
% the following are RD expenditures by Gov or Comp(anies)
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
SgSf=greenTotal./fossilTotal;

% => match Fpf/Y to shareProd: Af
% => outputFossil to match Emissions data: omegaa

pf =(pOil.*prodPetrol + pCoal.*prodCoal + pGas.*prodGas)./(prodPetrol + prodCoal + prodGas);
pfReal = pf./deflate*100; 

end