function [ForwardSim] = ForwardSim_GeyerModel_Explicit_AnkleOnly_vTwente(Data,Set,Reflex,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if ~isempty(varargin)
    RefCycle = varargin{1};
end

Bool_x0FromOpt = 0;
if length(varargin)>1
    Res = varargin{2};
    Bool_x0FromOpt = 1;
end

% get default settings
Set = GetDefaultSettings(Set);

tspan = Data.tspan;
fsSample = 1000;
tSample  = tspan(1):(1/fsSample):tspan(end);    % sampling time in optimization

DatSpline = ppval(Data.DataSpline,tSample)';
lMTdat = ppval(Data.lMT_Spline,tSample)';
dMdat = ppval(Data.dof(3).dM_Spline,tSample)';
BoolStance = DatSpline(:,strcmp(Set.headers,'BoolStance'));
BoolDS = DatSpline(:,strcmp(Set.headers,'DoubleSup_LeadRight'));
BoolSwing = DatSpline(:,strcmp(Set.headers,'BoolSwing'));

% filter moment data
iTa = find(strcmp(Set.headers,'Ta'));
iTk = find(strcmp(Set.headers,'Tk'));
iTh = find(strcmp(Set.headers,'Th'));
Tdat = DatSpline(:,[iTa iTk iTh]);
[a,b] = butter(4,6./(Set.MeshFreq*0.5),'low');
Tdat  = filtfilt(a,b,Tdat);
TSpline = spline(tSample,Tdat');

iGRFLz = strcmp(Set.headers,'GRFLz');
GRF_dat = DatSpline(:,iGRFLz);
GRF_dat = filtfilt(a,b,GRF_dat);
Fy_Norm  = GRF_dat./(Set.BodyMass*9.81);
Fy_Norm  = 0.5*tanh(10*(Fy_Norm-0.4))+0.5;

lMTdat = filtfilt(a,b,lMTdat);
lMTSpline = spline(tSample,lMTdat');

dMdat = filtfilt(a,b,dMdat);
dMSpline = spline(tSample,dMdat');

iStance = 1;
iDS = find(BoolDS>0.8,1,'first');
iSwing = find(BoolSwing>0.8,1,'first');
IndexVect = [iStance iDS iSwing];
tPhase = tSample(IndexVect);

[a,b] = butter(4,6./(fsSample*0.5),'low');
COMdSpline   = DatSpline(:,strcmp(Set.headers,'COMd'));
COMd         = filtfilt(a,b,COMdSpline);
COMSpline    = DatSpline(:,strcmp(Set.headers,'COM'));
COM          = filtfilt(a,b,COMSpline);


% Store input data in Misc structure (used in forward simulation)
Misc.Set            = Set;
Misc.Reflex         = Reflex; % reflex parameters
Misc.params         = Data.Muscle.params;
Misc.COM_TimeDelay  = Set.COM_TimeDelay;
Misc.ATendon        = [20 35]';
Misc.shift          = getShift(Misc.ATendon);

Misc.SP.lMT_Spline    = lMTSpline;
Misc.SP.dMS_qa_Spline = dMSpline;
Misc.SP.T_Spline      = TSpline;
Misc.SP.Fy_Spline     = spline(tSample,Fy_Norm);


%% Compute spline for deviation in COM kinematics

% deviation COM kinematics
if Set.PercGaitCycle
    tCOM = linspace(0,1,length(tSample));
    COMref = ppval(RefCycle.Spline_COM,tCOM);
    COMdref = ppval(RefCycle.Spline_COMd,tCOM);
    SP.COM      = COM - COMref';
    SP.COMd     = COMd - COMdref';
    if Set.COM_std
        [a,b] = butter(4,6./(fsSample*0.5),'low');
        % adapt input signal
        COMref_std = ppval(RefCycle.Spline_COMstd,tSample);
        COMdref_std = ppval(RefCycle.Spline_COMdstd,tSample);
        IndsBelow = abs(SP.COM(:,1))<COMref_std(:,1);
        SP.COM(IndsBelow,1) = 0;
        IndsBelow = abs(SP.COMd(:,1))<COMdref_std(:,1)*2;
        SP.COMd(IndsBelow,1) = 0;
        % filter input signal again to make it smooth
        SP.COMd     = filtfilt(a,b,SP.COMd);
        SP.COM      = filtfilt(a,b,SP.COM);
    elseif Set.ConstantSTD
        COMref_std = Set.COMref_std;
        COMdref_std = Set.COMdref_std;
        IndsBelow = abs(SP.COM(:,1))<COMref_std;
        SP.COM(IndsBelow,1) = 0;
        IndsBelow = abs(SP.COMd(:,1))<COMdref_std;
        SP.COMd(IndsBelow,1) = 0;
        % filter input signal again to make it smooth
        SP.COMd     = filtfilt(a,b,SP.COMd);
        SP.COM      = filtfilt(a,b,SP.COM);
    end
end
% compute time delayed COM feedback (% note: time delay is implemented in
% function for state derivative ExplicitDynamics_Geyer_AnkleOnly).
Misc.SP.COMdSpline = spline(tSample,SP.COMd);
Misc.SP.COMSpline = spline(tSample,SP.COM(:,1));



%% Integration implementation

if Bool_x0FromOpt
    % method to compute the initial state of the model
    % option 1: initial state based on result parameter estimation
    i0 = 1;
    x0 = [Res(1).a(:,i0);   Res(1).lMtilde(:,i0);   Res(1).Fpade(:,i0);     Res(1).lMpade(:,i0)]; % initial state
else
    i0 = 1;
    x0 = GetInitialState_Reflex_AnkleOnly_Fast(tSample(i0),Misc,Set);
end

% gradient decrease Force feedback during double support
Misc.BoolGradientDS = Set.BoolGradientDS;
Misc.G_gradient = Set.G_gradient;

% first integration step
tspan = [tPhase(1) tPhase(2)];

% pre-allocate output
tVect = nan(1000,1);
yVect = nan(1000,6);
PhaseVect = nan(1000,1);
ct = 0;
Ncycl = 1;
for i = 1:Ncycl
    % run integration - stance phase
    Misc.IndexPhase = 1; % stance phase =
    FunSel = @(t,y)ExplicitDynamics_Geyer_AnkleOnly(t,y,Misc);
    [sol] = ode23s(FunSel, tspan, x0);
    iSel = ct+1:ct+length(sol.x);
    tVect(iSel)     = sol.x';
    yVect(iSel,:)   = sol.y';
    PhaseVect(iSel) = Misc.IndexPhase;
    ct = ct+length(sol.x);
    
    % run integration - DS phase
    Misc.IndexPhase = 2; % stance phase
    x0 = sol.y(:,end);
    tspan = [tPhase(2+3*(i-1)) tPhase(3+3*(i-1))];
    Misc.t0_grad = tspan(1);
    FunSel = @(t,y)ExplicitDynamics_Geyer_AnkleOnly(t,y,Misc);
    [sol] = ode23s(FunSel, tspan, x0);
    
    iSel = ct+1:ct+length(sol.x);
    tVect(iSel)     = sol.x';
    yVect(iSel,:)   = sol.y';
    PhaseVect(iSel) = Misc.IndexPhase;
    ct = ct+length(sol.x);
    
    % run integration - swing phase
    Misc.IndexPhase = 0; % stance phase
    x0 = sol.y(:,end);
    if i <Ncycl
        tspan = [tPhase(3+3*(i-1)) tSwing(4+3*(i-1))];
    else
        tspan = [tPhase(3+3*(i-1)) tSample(end)];
    end
    FunSel = @(t,y)ExplicitDynamics_Geyer_AnkleOnly(t,y,Misc);
    [sol] = ode23s(FunSel, tspan, x0);
    iSel = ct+1:ct+length(sol.x);
    tVect(iSel)     = sol.x';
    yVect(iSel,:)   = sol.y';
    PhaseVect(iSel) = Misc.IndexPhase;
    ct = ct+length(sol.x);
    
    % update the intial state vector for next iteration
    if i<Ncycl
        x0 = sol.y(:,end);
        tspan = [tStance(i+1) tDS(i+1)];
    end
end
%
% remove pre-allocated
tVect(ct:end)     = [];
yVect(ct:end,:)   = [];
PhaseVect(ct:end) = [];

% extract the information
ForwardSim.t = tVect;
ForwardSim.phase = PhaseVect;
ForwardSim.a = yVect(:,1:2);
ForwardSim.lmTilde = yVect(:,3:4);
ForwardSim.Fpade = yVect(:,5);
ForwardSim.lMpade = yVect(:,6);


%% analyse the state

% selected muscles
iMuscle = [8 9];
ATendon = [20 35]';
shift   = getShift(ATendon);

% moment arms
dMs_qa   = ppval(Misc.SP.dMS_qa_Spline,ForwardSim.t)';    % moment arms ankle
lMT      = ppval(Misc.SP.lMT_Spline,ForwardSim.t)';        % Muscle tendon length

% muscle forces
ForwardSim.lMtilde_dot = zeros(length(ForwardSim.t),length(iMuscle));
ForwardSim.FT = zeros(length(ForwardSim.t),length(iMuscle));
for i =1:length(ForwardSim.t)
    [ForwardSim.lMtilde_dot(i,:),ForwardSim.FT(i,:)] = Explicit_ContractionDyn(ForwardSim.a(i,:)',...
        ForwardSim.lmTilde(i,:)',lMT(i,iMuscle)',Misc.params(:,iMuscle),ATendon,shift);
end

% joint moments
ForwardSim.Tankle = sum(ForwardSim.FT.*dMs_qa(:,iMuscle),2);
ForwardSim.Tid = ppval(Misc.SP.T_Spline,ForwardSim.t)';
ForwardSim.Tid_ankle = ForwardSim.Tid(:,1);
ForwardSim.Event = [];
ForwardSim.tPhase = tPhase;

end

