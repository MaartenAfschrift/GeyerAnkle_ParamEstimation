function [xdot] = ExplicitDynamics_Geyer_AnkleOnly(tSim,x,Misc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% select phase gait cycle
IndexPhase  = Misc.IndexPhase;

% unpack the vector with reflex gains
Set         = Misc.Set;
SP          = Misc.SP;
Sol         = Misc.Reflex.Sol;
Tib         = Misc.Reflex.Tib;
ATendon     = Misc.ATendon;
shift       = Misc.shift;
BoolGradientDS = Misc.BoolGradientDS;

%% Evaluate splines to get data

% Evaluate the Splines
lMT      = ppval(SP.lMT_Spline,tSim)';        % Muscle tendon length
Fy       = ppval(SP.Fy_Spline,tSim);
% if Set.COMfb
% COM_delay     = ppval(SP.COMSpline,tSim-Misc.COM_TimeDelay)';
COMd_delay    = ppval(SP.COMdSpline,tSim-Misc.COM_TimeDelay)';
% end

%% "unpack" state vector
nMuscles = 2;

% state
ct = 0;
a = x(ct+1:ct+nMuscles); ct = ct+nMuscles;
lMtilde = x(ct+1:ct+nMuscles); ct = ct+nMuscles;
yPade_Fb = x(ct+1:ct+1); ct = ct+1;
yPade_lM = x(ct+1:ct+1);

%% Explicit dynamics

iSol = 1;
iTib = 2;
iMuscle = [8 9];
iForce_FB = iSol; % indices adapted here because Fmus and lMtilde is smaller
ilM_FB = iTib; % indices adapted here because Fmus and lMtilde is smaller

% state derivative muscle dynamics
lMT = lMT(iMuscle)';
params  = Misc.params(:,iMuscle);
% [dlMtildedt,FT] = Explicit_ContractionDyn(a,lMtilde,lMT,params,ATendon,shift);
[dlMtildedt,FT] = Explicit_ContractionDyn_limitFMV(a,lMtilde,lMT,params,ATendon,shift);

% pade aproximation
Fmus_norm = FT./params(1,:)';
[Fmus_t,yPade_Fb_dot] = PadeApprox_order1_20ms(Fmus_norm(iForce_FB),yPade_Fb);
[lM_t,yPade_lM_dot] = PadeApprox_order1_20ms(lMtilde(ilM_FB),yPade_lM);

% compute muscle excitations
iSol_fMus = find(iSol == iForce_FB);  % index soleus in "-delayed" force
iTib_lM = find(iTib == ilM_FB); % index tibialis in "-delayed" fiber length

GRFConds = {'Geyer_Default_GRF','Geyer_COMd_GRF','Geyer_COMd_GRF_lim'};
if IndexPhase == 1 % stance phase
    % Soleus Reflex    
    if any(strcmp(Set.GainsEstimation,GRFConds))
        eSol = Sol.G.*Fy.*Fmus_t(iSol_fMus,:);
    else
        eSol = Sol.G.*Fmus_t(iSol_fMus,:);
    end
    COMfb = COMd_delay.*Sol.COMd;   
    eSol = eSol + COMfb;
    
    % Tibialis anterior reflex    
    eTib   = Tib.G.*(lM_t(iTib_lM)-Tib.loff) - Tib.GSol.*Fmus_t(iSol_fMus,:);
    COMfb = COMd_delay.* Tib.COMd;
    eTib = eTib + COMfb;
    
    
elseif IndexPhase == 2 % stance phase with pre-swing excitations double support
    % Soleus Reflex
    if BoolGradientDS
        tgrad = tSim-Misc.t0_grad;
        eSol = (Sol.G+Misc.G_gradient*tgrad).*Fmus_t(iSol_fMus,:);
    else
        if any(strcmp(Set.GainsEstimation,GRFConds))
            eSol = Sol.G.*Fy.*Fmus_t(iSol_fMus,:);    
        else
            eSol = Sol.G.*Fmus_t(iSol_fMus,:);    
        end
    end
    COMfb = COMd_delay.*Sol.COMd;
    eSol = eSol + COMfb;
    
    % Tibialis anterior reflex
    eTib   = Tib.G.*(lM_t(iTib_lM,:)-Tib.loff) - Tib.GSol.*Fmus_t(iSol_fMus,:);    
    COMfb = COMd_delay.* Tib.COMd;
    eTib = eTib + COMfb;
    
else % swing phase
    % Soleus Reflex
    eSol = 0;
    
    % Tibialis anterior reflex
    iTib_lM = find(iTib == ilM_FB); % index tibialis in "-delayed" fiber length
    eTib   = Tib.G.*(lM_t(iTib_lM,:)-Tib.loff);
end

% combine excitations in one vector
eFB = [eSol; eTib];

% allow saturation of values between 0 and 1
[eFB] = SaturateExcitation(eFB,15);

% add baseline activity after saturation
eBaseline = [Sol.e0 Tib.e0]';
eFB = eFB + eBaseline;

% muscle activation dynamics
dadt = ActivationDynamics(eFB,a,0.02,0.04,0.1);

% state derivative
xdot = [dadt; dlMtildedt; yPade_Fb_dot; yPade_lM_dot];


end

