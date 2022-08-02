function COMP = comp_CEV(sol1, sol2, varlist1, varlist2, symms, list, params, T, indic)
% this function compares policies! 

% input
% sols : cell of two policy results to be compared
%        the first one is the benchmark
%       i.e. the resulting CEV says: percentage change in consumption to go
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
disc=repmat(betaa, 1,T);
expp=0:T-1;
vec_discount= disc.^expp;

%- CEVs over 55 years
Wasterisk =  vec_discount*sol1(varlist1=='SWF',:)'; %benchmark swf

%- CEV over 55 periods 
if thetaa~=1
    U2_noCon  = vec_discount*(sol2(varlist2=='Utillab',:)'+sol2(varlist2=='Utilsci',:)');
    CEV   = transpose(ones(size(sol1(varlist1=='SWF',:))).*((Wasterisk + U2_noCon)/...
        (vec_discount*transpose(sol2(varlist2=='C',:)).^(1-thetaa))*(1-thetaa))^(1/(1-thetaa))-1);
    %- dynamic CEV
    CEVDy =transpose((sol1(varlist1=='SWF',:)+sol2(varlist2=='Utillab',:)+sol2(varlist2=='Utilsci',:)...
        ./(sol2(varlist2=='C',:).^(1-thetaa)).*(1-thetaa)).^(1/(1-thetaa))-1);
else
    SWF2  = vec_discount*sol2(varlist2=='SWF',:)';
    CEV   = transpose(ones(size(sol1(varlist1=='SWF',:))).*exp((Wasterisk-SWF2)/sum(vec_discount))-1);
    CEVDy = transpose(exp(sol1(varlist1=='SWF',:)-sol2(varlist2=='SWF',:))-1);
end

    
%- save results
COMP=eval(symms.comp); 

end