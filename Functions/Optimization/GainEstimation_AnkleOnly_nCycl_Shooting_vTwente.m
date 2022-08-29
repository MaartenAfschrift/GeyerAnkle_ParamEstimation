function [Solution,f_ForwardSim] = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(PhaseDat,Set,varargin)
%GainEstimation_AnkleOnly_nCycl_Shooting Optimize feedback gains of reflex driven model
% to track inverse dynamics joint moments using a shooting appraoch.
%   Input arguments:
%       (1) PhaseDat: Structure with MoCap information
%       (2) Set: Structure with all settings for the optimization.
%           (see GetDefaultSettings.m for details)

% Warning !
% This fast version assumes that you use an euler integration scheme with a
% step of 0.001s. (This is needed because the full integration scheme has been
% built before in casadi to save time (CreateCasadiFunctions_Shooting.m)

if ~isempty(varargin)
    RefCycle = varargin{1};
end


%% Default settings
Set = GetDefaultSettings(Set);

% double check if the sampling frequency is 1000Hz
if Set.MeshFreq ~= 1000
    warning('You can only use this function with an integration step of 1ms. Adapted Freq to 1000 Hz')
    Set.MeshFreq = 1000;
end

%% Evaluate splines at optimization mesh
%----------------------------------------

nDatSets = length(PhaseDat);
for np =1:nDatSets
        
    % time vector for integration
    tspan    = PhaseDat(np).tspan;
    tSim     = tspan(1):(1/Set.MeshFreq):tspan(end);    % sampling time in optimization
    
    % Evaluate the splines on the mesh points
    DatSpline = ppval(PhaseDat(np).DataSpline,tSim)';
    lMTSpline = ppval(PhaseDat(np).lMT_Spline,tSim)';
    vMTSpline = ppval(PhaseDat(np).vMT_Spline,tSim)';
    dMs_qa_Spline = ppval(PhaseDat(np).dof(3).dM_Spline,tSim)';
    
    iTa = find(strcmp(Set.headers,'Ta'));
    iTk = find(strcmp(Set.headers,'Tk'));
    iTh = find(strcmp(Set.headers,'Th'));
    TSpline = DatSpline(:,[iTa iTk iTh]);
    
    iqa = find(strcmp(Set.headers,'qa'));
    iqk = find(strcmp(Set.headers,'qk'));
    iqh = find(strcmp(Set.headers,'qh'));
    qSpline = DatSpline(:,[iqa iqk iqh]);
    
    iGRFLz = strcmp(Set.headers,'GRFLz');
    GRF_Spline = DatSpline(:,iGRFLz);
    
    iGRFLy = strcmp(Set.headers,'GRFLy');
    iGRFRy = strcmp(Set.headers,'GRFRy');
    COMdd_Spline2 = (DatSpline(:,iGRFLy) + DatSpline(:,iGRFRy))./Set.BodyMass;
    
    % add perturbation force
    Fpert = zeros(size(COMdd_Spline2));
    if ~any(isnan(PhaseDat(np).dtP)) && length(PhaseDat(np).dtP)>1 && ~(all(PhaseDat(np).dtP ==0))
        iPert = find(tSim >= PhaseDat(np).dtP(1) & tSim <= PhaseDat(np).dtP(2));
        Fpert(iPert) =PhaseDat(np).PMag;
    end
    %COMdd_Spline = COMdd_Spline2 + Fpert./Set.BodyMass;
    COMdd_Spline = Fpert./Set.BodyMass;
      
    % low pass filter input
    [a,b]   = butter(4,6./(Set.MeshFreq*0.5),'low');
    [a2,b2] = butter(4,20./(Set.MeshFreq*0.5),'low');
    SP.dMs_qa   = filtfilt(a,b,dMs_qa_Spline);
    SP.T        = filtfilt(a,b,TSpline);
    SP.lMT      = filtfilt(a,b,lMTSpline);
    SP.vMT      = filtfilt(a,b,vMTSpline);
    SP.q        = filtfilt(a,b,qSpline);
    SP.Fy       = filtfilt(a,b,GRF_Spline);
    SP.Fy_Norm  = SP.Fy./(Set.BodyMass*9.81);
    SP.Fy_Norm  = 0.5*tanh(10*(SP.Fy_Norm-0.4))+0.5;
    
    
    COMdSpline      = DatSpline(:,strcmp(Set.headers,'COMd'));
    SP.COMd         = filtfilt(a,b,COMdSpline);
    COMSpline       = DatSpline(:,strcmp(Set.headers,'COM'));
    COMSpline       = COMSpline - COMSpline(1);
    SP.COM          = filtfilt(a,b,COMSpline);
    SP.COMdd        = COMdd_Spline; % GRF are already filtered

    % Compute deviation in COM kinematics when using COM fb in the
    % controller
    if Set.PercGaitCycle % COM deviations from unperturbed reference trajectory
        
        % get average COM data at this datapoint
        tCOM = linspace(0,1,length(tSim));
        COMref = ppval(RefCycle.Spline_COM,tCOM);
        COMdref = ppval(RefCycle.Spline_COMd,tCOM);
        COMddref =  ppval(RefCycle.Spline_COMdd,tCOM);
        
        SP.COM_copy = SP.COM ;
        SP.COMd_copy = SP.COMd;
        SP.COMdd_copy = SP.COMdd;
        SP.COM      = SP.COM - COMref';
        SP.COMd     = SP.COMd - COMdref';
%         SP.COMdd    = SP.COMdd - COMddref';
        if Set.COM_std
            % adapt input signal
            COMref_std = ppval(RefCycle.Spline_COMstd,tCOM);
            COMdref_std = ppval(RefCycle.Spline_COMdstd,tCOM);
            IndsBelow = abs(SP.COM(:,1))<COMref_std';
            SP.COM(IndsBelow,1) = 0;
            IndsBelow = abs(SP.COMd(:,1))<COMdref_std'*2;
            SP.COMd(IndsBelow,1) = 0;
            % filter input signal again to make it smooth
            SP.COMd     = filtfilt(a,b,SP.COMd);
            SP.COM      = filtfilt(a,b,SP.COM);
        elseif Set.ConstantSTD
            COMref_std = Set.COMref_std;
            COMdref_std = Set.COMdref_std;
            COMddref_std = Set.COMddref_std;
            IndsBelow = abs(SP.COM(:,1))<COMref_std;
            SP.COM(IndsBelow,1) = 0;
            IndsBelow = abs(SP.COMd(:,1))<COMdref_std;
            SP.COMd(IndsBelow,1) = 0;
            IndsBelow = abs(SP.COMdd(:,1))<COMddref_std;
            SP.COMdd(IndsBelow,1) = 0;
            % filter input signal again to make it smooth
            SP.COMd     = filtfilt(a,b,SP.COMd);
            SP.COM      = filtfilt(a,b,SP.COM);
            SP.COMdd     = filtfilt(a2,b2,SP.COMdd);
        end
        
        if Set.BoolPlot
            BoolStance = DatSpline(:,strcmp(Set.headers,'BoolStance'));
            iStance = BoolStance>0.5;
            figure();
            subplot(2,3,1);
            plot(SP.COM_copy(iStance,1)); hold on; plot(COMref(iStance),'--k');
            subplot(2,3,2);
            plot(SP.COMd_copy(iStance,1)); hold on; plot(COMdref(iStance),'--k');
            subplot(2,3,3);
            plot(SP.COMdd_copy(iStance,1)); hold on; plot(COMddref(iStance),'--k');
            subplot(2,3,4);
            plot(SP.COM(iStance,1));
            xlabel('Frames'); ylabel('COM position error');
            subplot(2,3,5);
            plot(SP.COMd(iStance));
            xlabel('Frames'); ylabel('COM velocity error');
            subplot(2,3,6);
            plot(SP.COMdd(iStance));
            xlabel('Frames'); ylabel('COM acc error');
            suptitle('Deviation from COM trajectory used');
        end
    end
    % compute time delayed COM feedback
    % here we will use 50 ms
    COMdspline = spline(tSim,SP.COMd);
    tDelay = tSim - Set.COM_TimeDelay;
    SP.COMd_delay = ppval(COMdspline,tDelay);
    COMspline = spline(tSim,SP.COM(:,1));
    SP.COM_delay = ppval(COMspline,tDelay);
    COMddspline = spline(tSim,SP.COMdd);
    SP.COMdd_delay = ppval(COMddspline,tDelay);    
    
    BoolStance = DatSpline(:,strcmp(Set.headers,'BoolStance'));
    BoolDS = DatSpline(:,strcmp(Set.headers,'DoubleSup_LeadRight'));
    BoolSwing = DatSpline(:,strcmp(Set.headers,'BoolSwing'));
    
    iStance = 1;
    iDS = find(BoolDS>0.8,1,'first');
    iSwing = find(BoolSwing>0.8,1,'first');
    IndexVect = [iStance iDS iSwing];
    % append spline information to structure
    PhaseDat(np).SP         = SP;
    PhaseDat(np).tSim       = tSim;
    PhaseDat(np).IndexVect  = IndexVect;
    PhaseDat(np).nPhase     = 3;
    clear SP tSim IndexVect
end


%%  Reflex parameters depending on the settings
%----------------------------------------------

% general muscle information
iMuscle     = [8 9];

% feedback gains (depending on the parameter estimation case used)
ConstReflex = [];
if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
    OptReflex    = SX.sym('OptReflex',6);
elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
    OptReflex   = SX.sym('OptReflex',10);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
    OptReflex    = SX.sym('OptReflex',8);
elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
    OptReflex    = SX.sym('OptReflex',3);
    Sol.e0      = Set.Sol.e0;
    Tib.e0      = Set.Tib.e0;
    Tib.loff    = Set.Tib.loff;
    ConstReflex = [ Sol.e0  Tib.e0 Tib.loff];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
    OptReflex    = SX.sym('OptReflex',5);
    Sol.e0      = Set.Sol.e0;
    Tib.e0      = Set.Tib.e0;
    Tib.loff    = Set.Tib.loff;
    ConstReflex = [ Sol.e0  Tib.e0 Tib.loff];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
    OptReflex    = SX.sym('OptReflex',6);
    Sol.e0      = Set.Sol.e0;
    Tib.loff    = Set.Tib.loff;
    ConstReflex = [Sol.e0 Tib.loff];
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
    OptReflex    = SX.sym('OptReflex',7);
    Sol.e0       = Set.Sol.e0;
    ConstReflex = [Sol.e0];
elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
    OptReflex    = SX.sym('OptReflex',10);    
elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
    OptReflex    = SX.sym('OptReflex',12);
end

%%          Forward simulation
%--------------------------------------

% Note these functions are pre-built with the script: "CreateCasadiFunctions_Shooting.m"

% casadi libraries
import casadi.*
opti = casadi.Opti();

% initial state as optimization variables
x0 = SX.sym('x0',6,nDatSets);

% loop over the perturbations to pre-allocate matrices
Ntot = 0;
iDatSet = zeros(nDatSets,1);
for np =1:nDatSets
    nPhase = PhaseDat(np).nPhase;
    for  f =1:nPhase
        %   select mesh points
        i0 = PhaseDat(np).IndexVect(f);
        tSim = PhaseDat(np).tSim;
        if f<nPhase
            iend = PhaseDat(np).IndexVect(f+1);
        else
            iend = length(tSim);
        end
        % selected indexes
        iSel = i0:iend;
        N = length(iSel)-1;
        Ntot = Ntot + N -1;
    end
    if np>1
        Ntot = Ntot + 1;
    end
end
Ntot = Ntot + 1;

% build casadif functions using script: CrateCasadiFunctions_Shooting
% PathDefaultFunc = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\Geyer_InverseID\CasadiFuncs';
PathDefaultFunc = fullfile(Set.MainPath,'Functions\CasadiFuncs');
f_Int1000_Stance_Full = Function.load(fullfile(PathDefaultFunc,[Set.GainsEstimation 'f_Int1000_Stance_Full']));
f_Int1000_DS_Full = Function.load(fullfile(PathDefaultFunc,[Set.GainsEstimation 'f_Int1000_DS_Full']));
f_Int1000_Swing_Full = Function.load(fullfile(PathDefaultFunc,[Set.GainsEstimation 'f_Int1000_Swing_Full']));

% Ntot = Ntot + 1000;
% pre-allocate states, objective func and aux vars.
States              = SX(6,Ntot);
J                   = SX(Ntot,1);
TMusV               = SX(Ntot,1);
TidV                = nan(Ntot,1);
tV                  = nan(Ntot,1);
% loop over all datasets (nDatSets) and gait phases in a dataset(nPhase)
ctx = 1;
IndexNewSim = 1:nDatSets;
for np =1:nDatSets
    nPhase = PhaseDat(np).nPhase;
    iDatSet(np,1) = ctx;
    for  f =1:nPhase
        % set the initial state
        if f == 1
            if np >1
                ctx = ctx+1;
            end
            States(:,ctx) = x0(:,np);
            IndexNewSim(np) = ctx;
        end
        
        % index of gait cycle (1) stance, (2) DS pre-sing (0) swing
        % the data always start with a heelstrike (i.e. stance)
        IndexPhase = rem(f,3);
        tSim = PhaseDat(np).tSim;
        
        %   select mesh points
        i0 = PhaseDat(np).IndexVect(f);
        if f<nPhase
            iend = PhaseDat(np).IndexVect(f+1);
        else
            iend = length(tSim);
        end
        
        % selected indexes
        iSel = i0:iend;
        N = length(iSel)-1;
        
        % get spline input information
        lMT = PhaseDat(np).SP.lMT(iSel(1:end-1),iMuscle);
        COMd_delay = PhaseDat(np).SP.COMd_delay(iSel(1:end-1));
        dMs = PhaseDat(np).SP.dMs_qa(iSel(1:end-1),:);
        Tid = PhaseDat(np).SP.T(iSel(1:end-1),:);
        Fy = PhaseDat(np).SP.Fy_Norm(iSel(1:end-1));
        COM_delay = PhaseDat(np).SP.COM_delay(iSel(1:end-1));
        COMdd_delay = PhaseDat(np).SP.COMdd_delay(iSel(1:end-1));
        
        lMTin = [lMT; zeros(1000-N,2)];
        COMd_delayin = [COMd_delay'; zeros(1000-N,1)];
        COM_delayin = [COM_delay'; zeros(1000-N,1)];
        COMdd_delayin = [COMdd_delay'; zeros(1000-N,1)];
        dMsin = [dMs(:,8:9); zeros(1000-N,2)];
        Tidin = [Tid(:,1); zeros(1000-N,1)];
        Fyin = [Fy; zeros(1000-N,1)];
        
        if Set.ScaleTroque
            Tidin = Tidin./(Set.BodyMass*Set.Height).*(70*1.75);
        end
        
        % select the right casadif function to get the matrices with the
        % integration equations.
        if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || ...
                strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin');
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin');
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin',Fyin);
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin',Fyin);
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin',Fyin);
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin',Fyin);
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,ConstReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,States(:,ctx),lMTin',...
                    COM_delayin',COMd_delayin',COMdd_delayin',dMsin',Tidin',Fyin);
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,States(:,ctx),lMTin',...
                    COM_delayin',COMd_delayin',COMdd_delayin',dMsin',Tidin',Fyin);
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,States(:,ctx),lMTin',...
                    COMd_delayin',COMdd_delayin',dMsin',Tidin',Fyin);
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,States(:,ctx),lMTin',...
                    COMd_delayin',COMdd_delayin',dMsin',Tidin',Fyin);
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        else
            if IndexPhase == 1
                [StatesOut,Jout,TmusVOut] = f_Int1000_Stance_Full(OptReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin');
            elseif IndexPhase == 2
                [StatesOut,Jout,TmusVOut] = f_Int1000_DS_Full(OptReflex,States(:,ctx),lMTin',COMd_delayin',dMsin',Tidin');
            else
                [StatesOut,Jout,TmusVOut] = f_Int1000_Swing_Full(OptReflex,States(:,ctx),lMTin',dMsin',Tidin');
            end
        end
        States(:,ctx+1:ctx+N-1) = StatesOut(:,2:N);
        J(ctx+1:ctx+N-1) = Jout(2:N);
        TMusV(ctx+1:ctx+N-1) = TmusVOut(2:N);
        
        % aux outputs
        tV(ctx:ctx+N)       = tSim(iSel);
        TidV(ctx:ctx+N-1)   = Tid(:,1);
        
        ctx = ctx+N-1;
    end
end
iDatSet =[iDatSet; ctx];

% create casadi functions
f_ForwardSim = Function('f_ForwardSim',{OptReflex,x0},...
    {States,J,TMusV},{'OptReflex','x0'},{'States','J','TMusV'});
f_ForwardSim_J = Function('f_ForwardSim_J',{OptReflex,x0},...
    {J},{'OptReflex','x0'},{'J'});

% f_coll_map = f_ForwardSim.map(1,'thread',3);
if isfield(Set,'ParComp') && Set.ParComp
    if ~isfield(Set,'nCores')
        Set.nCores = 4;
    end
    f_coll_map_J = f_ForwardSim_J.map(1,'thread',Set.nCores);
else
    Set.ParComp =  false;
end


% initial guess reflex gains
if ~isempty(Set.ReflexGuess)
    IG_Reflexes    = Set.ReflexGuess;
else
    % frome the paper of Geyer (this is the matrix of the full lower limb model)
    IG_Reflexes = [   0.0122,     1.0566,     0.0042,    0.4658,    0.0010,    0.1912,...
        0.8134,    0.2145,    0.0060,    0.6539,   -0.0019,    0.1616,    0.0278,...
        0.9999,    0.6377,    0.0048,    0.0002,    0.0003,    0.0010,    0.0003,...
        0.7028,    0.1758,    0.0001,   10.8182,    0.0000,   10.8182,    0.1222,...
        0.0466,    0.0052,    0.0032,    0.0028];
    % depending on the settings
    if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        IG_Reflexes = IG_Reflexes([1:2 5:8]);
    elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
        IG_Reflexes = [IG_Reflexes([1:2 5:8]) zeros(1,4)];
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd')
        IG_Reflexes = [IG_Reflexes([1:2 5:8]) zeros(1,2)];
    elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
        IG_Reflexes = IG_Reflexes([1 6 7]);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
        IG_Reflexes = [IG_Reflexes([1 6 7]) 0 0];
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        IG_Reflexes = [IG_Reflexes([1 3 6 7]) 0 0];
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        IG_Reflexes = [0 IG_Reflexes([2 5:8])];
    else
        IG_Reflexes = zeros(size(OptReflex));
    end
end

% compute a reasonable initial state
x0_guess = zeros(np,6);
disp('Compute initial states');
for np =1:nDatSets
    i0 = 1;
    tSample                 = PhaseDat(np).tSim;
    Misc.SP.lMT_Spline      = PhaseDat(np).lMT_Spline;
    Misc.SP.dMS_qa_Spline   = PhaseDat(np).dof(3).dM_Spline;
    Misc.SP.T               = PhaseDat(np).SP.T(1,1);
    Misc.params             = PhaseDat(np).Muscle.params;
    x0                      = GetInitialState_Reflex_AnkleOnly_Fast(tSample(i0),Misc,Set);
    x0_guess(np,:)          = x0;
end
disp('End compute initial states');
disp('');



%%      Create optimization problem
%------------------------------------------

diary(Set.diaryName);
disp('Start parameter estimation:');


% intial state as optimization variable
if Set.OptInitState
    x0 = opti.variable(6,nDatSets);
else
    x0 = x0_guess';
end

% reflex parameters as optimimzation variables
if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
    OptReflexO = opti.variable(1,6);
elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
    OptReflexO = opti.variable(1,10);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
    OptReflexO = opti.variable(1,8);
elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
    OptReflexO = opti.variable(1,3);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
    OptReflexO = opti.variable(1,5);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
    OptReflexO = opti.variable(1,6);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
    OptReflexO = opti.variable(1,7);
elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
    OptReflexO = opti.variable(1,12);
elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
    OptReflexO = opti.variable(1,10);
end


% reflexes as optimization variables
%bounds on the reflexes
[Bounds] = GetBoundsReflexParam_AnkleOnly_vTwente_vShoot(Set,Set.e0_min,Set.G_min);
opti.subject_to(Bounds(:,1)' < OptReflexO< Bounds(:,2)');

% set the initial reflex
opti.set_initial(OptReflexO,IG_Reflexes);

% set the initial state
if Set.OptInitState
    % set intial guess
    opti.set_initial(x0,x0_guess');
    
    % based bounds initial state on the initial guess of x0
    Meanx0_guess = nanmean(x0_guess);
    x0_lowerb = Meanx0_guess*0.9;
    x0_upperb = Meanx0_guess*1.1;
    
    % some usefull limits on the lower and upper bounds
    x0_lowerLogical = [Set.e0_min Set.e0_min 0.4 0.4 0 0.02];
    iLowerLogical = x0_lowerb < x0_lowerLogical;
    x0_lowerb(iLowerLogical) = x0_lowerLogical(iLowerLogical);
    
    x0_upperLogical = [0.9 0.9 1.3 1.3 0.2 0.4];
    iUpperLogical = x0_upperb > x0_upperLogical;
    x0_upperb(iUpperLogical) = x0_upperLogical (iUpperLogical);
    
    % bounds on the initial state (we might want to changes this)
    for i = 1:1:nDatSets
        opti.subject_to(x0_lowerb' < x0(:,i) < x0_upperb');
    end
end

% set the initial reflex
opti.set_initial(OptReflexO,IG_Reflexes);

% expression for the objective function
if Set.ParComp
    J = f_coll_map_J(OptReflexO,x0);
else
    J = f_ForwardSim_J(OptReflexO,x0);
end

% minimize the feedback gains as well (regularisation)
JGeyer = sumsqr(OptReflexO).*0.00001;

% minimize muscle excitations ?
if Set.MinA
    [StatesSim,~,~] = f_ForwardSim(OptReflexO,x0);
    aSim = StatesSim(1:2,:);
    J_Act = sumsqr(aSim).*Set.MinA_Weight./length(aSim(1,:));
end

% minimize objective function
J_Av = sumsqr(J)./length(J);
Jtot = J_Av + JGeyer;
if Set.MinA
    Jtot = Jtot + J_Act;
end
opti.minimize(Jtot);

% test with initial guess
if Set.PlotGuess
    % run forward simulation with guess
    [Guess.x,Guess.J,Guess.Tmus] = f_ForwardSim(IG_Reflexes,x0_guess');
    % plot the IG solution
    hGuess = figure();
    subplot(2,3,1)
    plot(full(Guess.x(1:2,:)'));
    title('activation');
    subplot(2,3,2)
    plot(full(Guess.x(3:4,:)'));
    title('fiber length');
    subplot(2,3,3)
    plot(full(Guess.x(5:6,:)'));
    title('pade states');
    subplot(2,3,4:6)
    plot(full(Guess.Tmus),'b'); hold on;
    plot(TidV,'--k'); hold on;
end
% sovle the optimization problem
optionssol.ipopt.hessian_approximation  = 'limited-memory';
optionssol.ipopt.nlp_scaling_method     = 'gradient-based';
optionssol.ipopt.linear_solver          = 'mumps';
optionssol.ipopt.tol                    = 1e-4;
optionssol.ipopt.max_iter               = 2000;
% optionssol.jit = true;
% optionssol.compiler = 'shell';
% optionssol.jit_options.verbose = true;
% optionssol.jit_options.flags = {'/O2'};
[w_opt, stats,lbx,ubx,llb,uub] = solve_NLPSOL(opti,optionssol);

% extract the optimal solution
[Solution.ReflexOpt, Solution.x0] = UnpackSolutionShooting_AnkleOnly(w_opt,Set,nDatSets);
if ~Set.OptInitState
    Solution.x0 = x0;
end

if Set.PlotSolution
    % run forward simulation with guess
    [Guess.x,Guess.J,Guess.Tmus] = f_ForwardSim(Solution.ReflexOpt,Solution.x0);
    % plot the IG solution
    if Set.PlotGuess  == true
        figure(hGuess);
    else
        figure();
    end
    subplot(2,3,1); hold on;
    plot(full(Guess.x(1:2,:)'),'r');
    title('activation');
    subplot(2,3,2); hold on;
    plot(full(Guess.x(3:4,:)'),'r');
    title('fiber length');
    subplot(2,3,3); hold on;
    plot(full(Guess.x(5:6,:)'),'r');
    title('pade states');
    subplot(2,3,4:6); hold on;
    plot(full(Guess.Tmus),'r'); hold on;
    plot(TidV,'--k'); hold on;
end

%% Simulate with optimal solution
%--------------------------------

[xSol,JSol,TmusSol] = f_ForwardSim(Solution.ReflexOpt,Solution.x0);

% unpack solution
Solution.x = full(xSol);
Solution.J = full(JSol);
Solution.Tmus = full(TmusSol);

% also append the ID moments
Solution.Tid = TidV;
Solution.t = tV;
Solution.iDatSet = iDatSet;

% Create structure for Reflex gains
[Solution.Reflex] =StructureReflexGains_Shooting(Set,Solution.ReflexOpt,ConstReflex);


%% Identification analysis

% Parameter identification analysis through covariance matrix of jacobian
% of the cost function with respect to the parameters.
if Set.Vlugt2010
    % function for the objective
    J = f_ForwardSim_J(OptReflex,Solution.x0);
    % create the jacobian
    Jac = jacobian(J,OptReflex);
    f_Jac = Function('f_Jac',{OptReflex},{Jac});
    
    % evaluate jacobian and objective numerically
    Jsim = full(f_ForwardSim_J(Solution.ReflexOpt,Solution.x0));
    JacSim = full(f_Jac(Solution.ReflexOpt));
    
    N = length(J);
    Vlugt.P = (1/N) * inv(((JacSim.') * JacSim)) * (Jsim' * Jsim);
    np = length(Vlugt.P);
    Vlugt.P_vis = nan(np,np);
    for i = 1:np
        for b = 1:np
            Vlugt.P_vis(i,b) = Vlugt.P(i,b)/ sqrt(Vlugt.P(i,i)*Vlugt.P(b,b));
        end
    end
    if Set.PlotVlugt
        PlotVlugt(Vlugt.P_vis,Set.GainsEstimation)
    end
    Vlugt.GainID = Set.GainsEstimation;
    Solution.Vlugt = Vlugt;
end

%% Store the settings and close diary
Solution.Set = Set;
Solution.stats = stats;
for np =1:nDatSets
    Solution.Phase(np).IndexVect = PhaseDat(np).IndexVect;
end

% end the log
diary off

end

