function [x0,x0_dot,R] = GetInitialState_Reflex_AnkleOnly_Fast(tSim,Misc,Set)
%GetInitialState_Reflex_AnkleOnly_Fast finds initial states that satisfies dynamics
% with minimal difference between muscle activation and excitation from
% feedback control

% main concept:
% find initial state with state deriative equals 0 and minimizes objective:
% J = (eFb - a)^2 + w1*(Tid-Tmus)^2 + w2*a^2


%% Evaluate splines to get data

ConvToRad = 1;
if ~isfield(Set,'InDeg') && Set.InDeg == 1
    ConvToRad = pi./180;
end

% evaluate the splines
dMs_qa   = ppval(Misc.SP.dMS_qa_Spline,tSim)';  % moment arms ankle
lMT      = ppval(Misc.SP.lMT_Spline,tSim)';        % Muscle tendon length
if isfield(Misc.SP,'T_Spline')
    T        = ppval(Misc.SP.T_Spline,tSim)';           % joint moments
else
    T = Misc.SP.T;
end

%% optimization variables
import casadi.*

opti = casadi.Opti();

% states
a = opti.variable(2,1);
lMtilde = opti.variable(2,1);
opti.set_initial(a,0);
opti.set_initial(lMtilde,1);

opti.subject_to(0.03 < a < 1);
opti.subject_to(0.2 < lMtilde < 1.8);

% static situation
lMtilde_dot = zeros(2,1);

%% Select only the soleus and gastrocnemius muscle in this model 
iMuscle = [8 9];

% projected fiber length
lMo = Misc.params(2,iMuscle);
alphao = Misc.params(4,iMuscle);
w = (lMo.*sin(alphao))';
lM = lMtilde.*lMo';
lM_proj = sqrt(lM.^2 - w.^2);

% force equilibrium
ATendon = [ones(1,6)*35 20 20 35]';
shift   = getShift(ATendon);
[c_HillDiff,Fmus] = ForceEquilibrium_lMtildeState(a,lMtilde,lMtilde_dot,...
    lM_proj,lMT(iMuscle)',Misc.params(:,iMuscle)',ATendon(iMuscle),shift(iMuscle));

% constraints
opti.subject_to(c_HillDiff == 0);

% compute the joint moments
TMus = sum(dMs_qa(iMuscle)'.*Fmus);

% pade approximation based on dFpade/dt = 0
Fiso = Misc.params(1,:);
Fmus_norm = Fmus(1,:)./Fiso(8)';
Fpade = 16*Fmus_norm./100;
LMpade = 16*lMtilde(2)./100;
Fmus_t = PadeApprox_order1_20ms_implicit(Fmus_norm,Fpade,0);
lM_t = PadeApprox_order1_20ms_implicit(lMtilde(2),LMpade,0);

%% objective function
J = sumsqr(T(1)-TMus);
J = J + 0.001*sumsqr(Fmus_t-Fmus(1));
J = J + 0.001*sumsqr(lM_t-lMtilde(2));
J = J + 0.1*sumsqr(a);
opti.minimize(J);

% solve optimization problem
optionssol.ipopt.nlp_scaling_method     = 'gradient-based';
optionssol.ipopt.linear_solver          = 'mumps'; % you might have to use mumps if ma57 is not installed
optionssol.ipopt.tol                    = 1e-5;
optionssol.ipopt.max_iter               = 10000;
optionssol.ipopt.print_level            = 0;
opti.solver('ipopt',optionssol);
sol = opti.solve();

% get optimal values
R.a = value(sol,a);
R.lMtilde = value(sol,lMtilde);
R.Fmus = value(sol,Fmus); 

% get the pade results
R.Fpade = value(sol,Fpade);
R.LMpade = value(sol,LMpade);

% initial state
x0 = [R.a; R.lMtilde; R.Fpade; R.LMpade];
x0_dot = zeros(22,1);

end

