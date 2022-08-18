function [eout] = SaturateExcitation(e,b)
%SaturateExcitation smooth saturatin of excitations between 0 and 1



e1      = e.*(tanh(b*e)*0.5+0.5);
e2      = (e1-1).*(tanh(b*(e1-1))*0.5+0.5);
eout    = e1-e2;


end

