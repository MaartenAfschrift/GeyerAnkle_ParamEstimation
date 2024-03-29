% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

function [err, FT, Fpe, FMltilde, FMvtilde, cos_alpha] = ForceEquilibrium_lMtildeState(a,lMtilde,vMtilde,lM_projected,lMT,params,kT,shift)

FMo = params(:,1);
lMo = params(:,2);
lTs = params(:,3);
vMtildemax = params(:,5);

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
lT = lMT - lM_projected;
lTtilde = lT./lTs;

% Tendon force-length characteristic
fse = (exp(kT.*(lTtilde - 0.995)))/5-0.25+shift;

% get muscle force-length characteristic
% lMtildeTemp = (lMtilde-1).*0.75+1;
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,vMtildemax);
% Fpe = Fpe.*0.1;
% Active muscle force
d = 0.00001; % damping coefficient
Fce = a.*FMltilde.*FMvtilde +d*vMtilde;

% Muscle force
FM = Fce+Fpe;
% Tendon force
FT = FMo.*fse;

% Equilibrium between muscle and tendon forces
% Fm*cos(alpha) = Ft
cos_alpha = (lMT-lT)./lM;
err =  FM.*cos_alpha-fse;

end