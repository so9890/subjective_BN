function allo_trans=trans_allo_out(indexx, solution, params, listt, indic)

% function to transform solution of fsolve or vpasolve with transformed
% model.

% input
% indexx: structure of logical vector ==1 if respective transformation has
%           to be performed
% solution: raw, untransformed solution from fsolve or vpasolve
% params: structure or double of parameter values: may be numeric

list.params=listt;
read_in_params;
%-- transformation 

allo_trans=solution;
allo_trans(indexx.sqr)=solution(indexx.sqr).^2;
allo_trans(indexx.exp)=exp(solution(indexx.exp));
allo_trans(indexx.lab)=params(listt=='upbarH')*indic.minn./(1+exp(solution(indexx.lab)));
if isfield(indexx, 'BN')
    allo_trans(indexx.BN)=B./(1+exp(solution(indexx.BN)));
end
if isfield(indexx, 'BNh')
    allo_trans(indexx.BNh)=Bh./(1+exp(solution(indexx.BNh)));
end
if isfield(indexx, 'BNl')
    allo_trans(indexx.BNl)=Bl./(1+exp(solution(indexx.BNl)));
end
allo_trans(indexx.oneab)=1./(1+exp(solution(indexx.oneab)));

end
