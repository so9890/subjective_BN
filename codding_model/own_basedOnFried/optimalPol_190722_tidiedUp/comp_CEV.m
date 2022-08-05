function COMP = comp_CEVv(sol1, sol2, varlist1, varlist2, symms, list, params, T, indic)
% this function compares policies! 

% input
% sols : cell of two policy results to be compared
%        the first one is the benchmark
%       i.e. the resulting CEVv says: percentage change in consumption to go
%       from sol2 to sol1

read_in_params;

% consumption equivalence

% 
% for count=1:2
%     varrs = sols{count}; 
%     liss= varrlist{count}; % pick adequate variable list
%     
%     C=varrs(liss=='C',:)';
%     hh = varrs(liss=='hh',:)';
%     hl = varrs(liss=='hl',:)';
%     sg = varrs(liss=='sg',:)';
%     sff = varrs(liss=='sff',:)';
%     sn = varrs(liss=='sn',:)';
%     
%     %- utility functions
% 
%     if thetaa~=1
%         Utilcon = (C.^(1-thetaa))./(1-thetaa);
%     elseif thetaa==1
%         Utilcon = log(C);
%     end
% 
%     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
% 
%     if indic.sep==0
%           Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
%      else
%           Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
%     end
%     
%     %- save stuff
%     
% end

%- create discount vector
%a) for version with PV: have to include last optimization period
disc=repmat(betaa, 1,T);
expp=0:T-1;
vec_discountPV= disc.^expp;

%b) for version without PV (exclude last period)
disc=repmat(betaa, 1,T-1);
expp=0:T-2;
vec_discountnoPV= disc.^expp;

%- CEVvs over 55 years
Wasterisk_noPV =  vec_discountnoPV*sol1(varlist1=='SWF',1:end-1)'; %benchmark swf
Wasterisk_withPV =  vec_discountPV*sol1(varlist1=='SWF',:)'; %benchmark swf
gammay = sol2(varlist2=='Y',end)/sol2(varlist2=='Y',end-1)-1; % approximate growth rate

%- CEVv over 55 periods 
if thetaa~=1
    U2_noConnoPV  = vec_discountnoPV*(sol2(varlist2=='Utillab',1:end-1)'+sol2(varlist2=='Utilsci',1:end-1)');
    U2_noConwithPV  = vec_discountPV*(sol2(varlist2=='Utillab',:)'+sol2(varlist2=='Utilsci',:)');
    PV2_noCon = betaa^T/(1-betaa).*(sol2(varlist2=='Utillab',end)'+sol2(varlist2=='Utilsci',end)'); 
    PV2_Con = betaa^T/(1-betaa*(1+gammay)^(1-thetaa))*sol2(varlist2=='C',end)^(1-thetaa)/(1-thetaa); 
    
    CEVv   = transpose(ones(size(sol1(varlist1=='SWF',:)))).*(((Wasterisk_noPV + U2_noConnoPV)/...
        (vec_discountnoPV*transpose(sol2(varlist2=='C',1:end-1)).^(1-thetaa))*(1-thetaa))^(1/(1-thetaa))-1);
    
    CEVvPV =  transpose(ones(size(sol1(varlist1=='SWF',:)))).*(((Wasterisk_withPV + sol1(varlist1=='PV',1)+ U2_noConwithPV+PV2_noCon)/...
        (vec_discountPV*transpose(sol2(varlist2=='C',:)).^(1-thetaa)/(1-thetaa)+PV2_Con))^(1/(1-thetaa))-1);
    %- dynamic CEVv
    CEVvDy =transpose((sol1(varlist1=='SWF',:)+sol2(varlist2=='Utillab',:)+sol2(varlist2=='Utilsci',:)...
        ./(sol2(varlist2=='C',:).^(1-thetaa)).*(1-thetaa)).^(1/(1-thetaa))-1);
else
    SWF2_noPV  = vec_discountnoPV*sol2(varlist2=='SWF',1:end-1)';
    SWF2_withPV  = vec_discountPV*sol2(varlist2=='SWF',:)';
    CEVv   = transpose(ones(size(sol1(varlist1=='SWF',:)))).*(exp((Wasterisk_noPV-SWF2_noPV)/sum(vec_discountnoPV))-1);
    CEVvPV   = transpose(ones(size(sol1(varlist1=='SWF',:)))).*(exp((Wasterisk_withPV+ sol1(varlist1=='PV',1)-SWF2_withPV- sol2(varlist2=='PV',1))...
        /(sum(vec_discountPV)+betaa^T./(1-betaa.*(1+gammay)^(1-thetaa))))-1);
    CEVvDy = transpose(exp(sol1(varlist1=='SWF',:)-sol2(varlist2=='SWF',:))-1);
end

    
%- save results
COMP=eval(symms.comp); 

end