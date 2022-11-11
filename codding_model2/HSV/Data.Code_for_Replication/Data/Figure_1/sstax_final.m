%*********************************
%PROGRAM SSTAX
%**********************************
clc;
clear all; 

%% Regression coefficients
hbetamat = [ 7.065382 .1091357  -.0004842  -8.98e-06;
             7.723434 .1478196  -.0025077  .0000133;
             3.042195 .4816315 -.0093977 .0000592]; 

wbetamat = [ 10.60823 -.2172776 .0081305  -.0000826;
             7.532779  .1296373 -.0022193 .0000122;
             7.185836  .1909089  -.0036502  .0000224];

%% Age profiles
J0 = 25;
Jmax = 60;
Jret = 65;
J    = 100;
hy = zeros(3,Jmax-J0+1);
wy = zeros(3,Jmax-J0+1);
for edu = 1:3
         for j = J0:Jmax
             age(j-(J0-1))  = j;
            hy(edu,j-(J0-1)) = exp(hbetamat(edu,1) + hbetamat(edu,2)*j + hbetamat(edu,3)*j^2 + hbetamat(edu,4)*j^3);
           wy(edu,j-(J0-1)) = exp(wbetamat(edu,1) + wbetamat(edu,2)*j + wbetamat(edu,3)*j^2 + wbetamat(edu,4)*j^3);
         end
end

%% Compute AIME
hy_aime    = zeros(3,Jmax-(J0-1));
hy_aimelag = zeros(3,Jmax-(J0-1));
hy_aimejk  = zeros(3,1);
wy_aime    = zeros(3,Jmax-(J0-1));
wy_aimelag = zeros(3,Jmax-(J0-1));
wy_aimejk  = zeros(3,1);
for edu = 1:3
         hy_aime(edu,1)    = hy(edu,1)/ (Jmax-(J0-1));
         hy_aimelag(edu,1) = 0;
         wy_aime(edu,1)    = wy(edu,1)/ (Jmax-(J0-1));
         wy_aimelag(edu,1) = 0;
         for j = J0+1:Jmax      
            hy_aime(edu,j-(J0-1)) = hy_aime(edu,j-1-(J0-1)) + hy(edu,j-(J0-1))/ (Jmax-(J0-1));
            hy_aimelag(edu,j-(J0-1)) = hy_aime(edu,j-1-(J0-1));
            wy_aime(edu,j-(J0-1)) = wy_aime(edu,j-1-(J0-1)) + wy(edu,j-(J0-1))/ (Jmax-(J0-1));
            wy_aimelag(edu,j-(J0-1)) = wy_aime(edu,j-1-(J0-1));
         end
         hy_aimejk(edu)  = hy_aime(edu,Jmax-(J0-1));
         wy_aimejk(edu)  = wy_aime(edu,Jmax-(J0-1));
end
hy_aime    = hy_aime/12;
hy_aimelag = hy_aimelag/12;
hy_aimejk  = hy_aimejk/12;
wy_aime    = wy_aime/12;
wy_aimelag = wy_aimelag/12;
wy_aimejk  = wy_aimejk/12;


%% Compute demographic factor 
% Life tables https://www.ssa.gov/oact/STATS/table4c6.html
data  = xlsread('lifetables.xlsx','a2:j42');
hsurvrate = (1 - data(:,2));
wsurvrate = (1 - data(:,5));
R     = 1.04;

%DPV from age Jret onward
hfactor = 0;
wfactor = 0;
for j = Jret:J
     hfactor = hfactor + hsurvrate(j-Jmax+1)/R^(j-Jret);
     wfactor = wfactor + wsurvrate(j-Jmax+1)/R^(j-Jret);
end
hadjfactor = zeros(Jmax-(J0-1),1);
wadjfactor = zeros(Jmax-(J0-1),1);
for j = J0:Jmax
hadjfactor(j-(J0-1),1) = hfactor*1/R^(Jret-j);
wadjfactor(j-(J0-1),1) = wfactor*1/R^(Jret-j);
end

%% Reformat for STATA
hy_aimevec = reshape(hy_aime',[3*(Jmax-(J0-1)),1]);
hy_aimelagvec = reshape(hy_aimelag',[3*(Jmax-(J0-1)),1]);
wy_aimevec = reshape(wy_aime',[3*(Jmax-(J0-1)),1]);
wy_aimelagvec = reshape(wy_aimelag',[3*(Jmax-(J0-1)),1]);

hy_aimejkvec = kron(hy_aimejk, ones(size(age)));
hy_aimejkvec = reshape(hy_aimejkvec',[3*(Jmax-(J0-1)),1]);
wy_aimejkvec = kron(wy_aimejk, ones(size(age)));
wy_aimejkvec = reshape(wy_aimejkvec',[3*(Jmax-(J0-1)),1]);

agevec = [age'; age'; age'];
eduvec = kron([ 1; 2; 3], ones(size(age)));
eduvec = reshape(eduvec',[3*(Jmax-(J0-1)),1]);

hadjfactorvec = [hadjfactor;hadjfactor;hadjfactor];
wadjfactorvec = [wadjfactor;wadjfactor;wadjfactor];

%This matrix can be merged with the stata dataset and has all info we need
hdatamatlab= [eduvec agevec hy_aimevec hy_aimelagvec hadjfactorvec hy_aimejkvec];
wdatamatlab= [eduvec agevec wy_aimevec wy_aimelagvec wadjfactorvec wy_aimejkvec];

%Now, print datamatlab in matlab, copy and paste into the excel file with
%same name, and from there copy into statafile. All manual :-(




