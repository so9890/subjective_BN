function [c, ceq]=sym_fmincon(x, indic)

c=[];
if indic.target==0
    ceq=Ram_Model_notarget_1905(x);
else
    ceq=Ram_Model_target_1905(x);
end
end
