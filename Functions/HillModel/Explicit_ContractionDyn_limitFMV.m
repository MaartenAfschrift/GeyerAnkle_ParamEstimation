function [dlMtildedt,FT] = Explicit_ContractionDyn_limitFMV(a,lMtilde,lMT,params,ATendon,shift,varargin)
% Explicit_ContractionDyn_limitFMV Computes time derivative of the normalised muscle
% fiber length. This is an adapted version that limits the maximal fiber velocity
% inputs:
%   a: muscle activation
%   lMtilde: normalised muscle fiber length
%   lMT: muscle-tendon length
%   params: muscle parameters
%   ATendon: value for tendon compliance
%   shift: shift in tendon F/L curve

% Create vector with parameters
FMo = params(1,:)';
lMo = params(2,:)';
lTs = params(3,:)';
alphao = params(4,:)';
vMtildemax = params(5,:)';

% Hill model geometric relationship
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
lT = lMT - sqrt((lM.^2 - w.^2));
lTtilde = lT./lTs;

% compute the tendon force
fse = (exp(ATendon.*(lTtilde - 0.995)))./5-0.25+shift;

% force-length properties
Faparam = [0.814483478343008;1.05503342897057;0.162384573599574;0.0633034484654646; ...
    0.433004984392647;0.716775413397760;-0.0299471169706956;0.200356847296188];

b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilde = FMtilde1+FMtilde2+FMtilde3;

% passive force
Fpparam = [-0.995172050006169; 53.598150033144236];
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 0.10e1) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

% derive F/v location
FMce = fse.* lM ./(lMT-lT) - Fpe;
FMvtilde = FMce./(a.*FMltilde);

% limit FMvtilde
off = 1.4; % FMvtilde tops of at this value
b = 100; % smooth transition
FMvtilde     = FMvtilde - (FMvtilde-off ).*(tanh(b*(FMvtilde-off ))*0.5+0.5);

% force-velocity properties
Fvparam = [-0.318323436899128   -8.149156043475250 ,...
    -0.374121508647860    0.885644059915004];

e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

vMtilde_Curve = 1/e2*(sinh((FMvtilde-e4)/e1)-e3);
vM = vMtilde_Curve .*vMtildemax .* lMo;
dlMtildedt = vM ./ lMo;   %  = vMtilde;
FT = fse.*FMo;

end

