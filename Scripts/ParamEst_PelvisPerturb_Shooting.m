
%% Gain estimation data Twente
%-----------------------------------

clear all; close all; clc;

% path information
Set.MainPath        = 'C:\Users\Maarten\Documents\Software\Sim\GeyerAnkle_ParamEstimation'; % Path to main folder of this repo

% relative paths
Set.PolyPath        = fullfile(Set.MainPath,'OsimModel\MuscleAnalysis');
Set.ModelPath       = fullfile(Set.MainPath,'OsimModel\Models','gait_2D_Strength.osim');
Set.ResultsFolder   = fullfile(Set.MainPath,'Results','ResParamID');
Set.datapath        = fullfile(Set.MainPath,'Data');    % path with data of Vlutters 2018 organized in format for parameter estimation

% add paths
addpath(genpath(fullfile(Set.MainPath,'Functions')));
import casadi.*

% Select Dataset
Set.TwenteData = 1;   % Boolean for selecting data Twente

% subject IDs5
% Set.SubjIDs = 5; 
Set.SubjIDs = [4 5 7 8 10];

% Body mass
Set.BodyMass = 70;          % in kg

% angles in radians or deg to evaluate polynomials
Set.BoolRad = false;

% set info from the model
Set.MuscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};

% initial guess based on previous solution
Set.PathIG = fullfile(pwd,'Reflex_AddHip.mat');

% boolean to plot figures during optimization
Set.BoolPlot = 0;

% scale tohe optimizaton problem ?
Set.ScaleGains = 1;

% minimal e0 (reflex gain)
Set.e0_min = 0.02;

% Gradient change force feedback between heelstrike contra-lateral leg and
% toe-off event (i.e. during double support)
Set.BoolGradientDS = false;
Set.G_gradient = -6;

% offset in ankle angle
Set.qa_offset = 0.04; % in radians)

% optimizize initial state as well ?
Set.OptInitState = false;

% Vlugt (covariance analysis of optimization parameters)
Set.Vlugt2010 = false;
Set.PlotVlugt = false;

% constant "band around COM data" to detect deviation in COM movement
Set.ConstantSTD = true; % based on STD unperturbed walking
Set.COMref_std = 0.04; % 2 cm deviation from average
Set.COMdref_std = 0.06; % 6 cm/s deviation in COM velocity from average
Set.COMddref_std = 0.01; % 1cm/s2 deviation in COM acceleration from average

% perturbation settings
Set.COMfb = 1;
Set.COM_std = false;
Set.PercGaitCycle = true;
Set.COM_TimeDelay = 0.06; % time delay COM feedback [s]

% create results directory
if ~isfolder(Set.ResultsFolder)
    mkdir(Set.ResultsFolder);
end

% information subjects Twente Dataset
% mass and height info subjects
Mass    = [ 66 89 53 62 56 80 76.3 56.6 57 66.6];
Height  = [1.79 1.92 1.62 1.82 1.80 1.94 1.89 1.73 1.78 1.70];

% scale torque to reference person
Set.ScaleTroque = true;

% download the dataset if needed
if ~isfolder(Set.datapath)
    Set.datapath = GetDataVlutters2018();
end

%% Create a Batch for simulations
% 

% Batch with all the results of the paper
% Slow = slow walking (0.62 m/s) and fast is normal walking speed (1.25m/s)
Batch.Speeds ={'Slow','Slow'};
% 1->4: backward increasing magnitude   5->8: forward increasing magnitude
Batch.iPertSel={1:8,...
    [1 2 4 5 6 8]};
Batch.nUnpSel = {8,8}; % number of unperturbed gait cycles
Batch.nP = {2,2};      % number of perturbed gait cycles

Batch.OutNames = {'Slow_AllDir_AllMag',...
    'Slow_Valid_Excl37'};


% % --- example of a batch with a hughe set of analysis
% % Slow = slow walking (0.62 m/s) and fast is normal walking speed (1.25m/s)
% Batch.Speeds ={'Slow','Slow','Slow','Fast','Fast','Fast',...
%     'Slow','Slow','Slow','Slow',...
%     'Fast','Fast','Fast','Fast',...
%     'Slow','Fast'};
% % 1->4: backward increasing magnitude   5->8: forward increasing magnitude
% Batch.iPertSel={1:8,1:4,5:8,1:8,1:4,5:8,...
%     1:2,3:4,6:7,7:8,...
%     1:2,3:4,6:7,7:8,...
%     [1 2 4 5 6 8],[1 2 4 5 6 8]};
% Batch.nUnpSel = {8,8,8,8,8,8,...
%     4,4,4,4,...
%     4,4,4,4,...
%     8,8}; % number of unperturbed gait cycles
% Batch.nP = {2,2,2,2,2,2,...
%     3,3,3,3,...
%     3,3,3,3,...
%     2,2};      % number of perturbed gait cycles
% 
% Batch.OutNames = {'Slow_AllDir_AllMag','Slow_Backward_AllMag','Slow_Forward_AllMag',...
%     'Fast_AllDir_AllMag','Fast_Backward_AllMag','Fast_Forward_AllMag',...
%     'Slow_AllDir_Mag12','Slow_AllDir_Mag34','Slow_AllDir_Mag56','Slow_AllDir_Mag78',...
%     'Fast_AllDir_Mag12','Fast_AllDir_Mag34','Fast_AllDir_Mag56','Fast_AllDir_Mag78',...
%     'Slow_Valid_Excl37','Fast_Valid_Excl37'};

% % -- minimal example
% Batch.Speeds ={'Slow'};
% %1->4: backward increasing magnitude   5->8: forward increasing magnitude
% Batch.iPertSel={1:8};
% % number of unperturbed gait cycles
% Batch.nUnpSel = {8};
% % number of repititions of unique perturbations
% Batch.nP = {2};
% % name output file
% Batch.OutNames = {'Slow_AllDir_AllMag'};


% number of batches
Batch.N = length(Batch.Speeds);                   % number of batches
Batch.SubjError = zeros(10,Batch.N);
SetO = Set; % copy of original settings

% check if the casadifunctions are compiled 
if ~isfolder(fullfile(Set.MainPath,'Functions','CasadiFuncs')) || ...
        ~exist(fullfile(Set.MainPath,'Functions','CasadiFuncs','GeyerDefaultf_Int1000_Swing_J'),'file')
    error(['You have to run the script ' fullfile(Set.MainPath,'Scripts','ParamEst_PelvisPush_Shooting') , ...
        'first to create the expression graphs with the system dynamics']);
end

for ib = 1:Batch.N

    % use original settings
    Set = SetO;

    % Select Data
    Set.WalkSpeed = Batch.Speeds{ib}; % walking speed
    Set.iPertSel = Batch.iPertSel{ib};
    Set.nUnp = Batch.nUnpSel{ib};
    Set.nP = Batch.nP{ib};

    % General OutputName
    Set.OutName_Gen = Batch.OutNames{ib};

    %% Loop over alls subjects
    for s = Set.SubjIDs
        try
        % copy from the original settings
        Set.sID = s;
        Set.BodyMass = Mass(s);          % in kg
        Set.Height = Height(s);          % in kgHeight

        if s == 1 || s == 6
            Set.qa_offset = 0.1; %larger offset in these subjects. problem with callibration
        else
            Set.qa_offset = SetO.qa_offset; %larger offset in these subjects. problem with callibration
        end

        %% Load the data of this subject
        Set.SubjName = ['pp_' num2str(Set.sID)];
        load(fullfile(Set.datapath,Set.SubjName,'RefWalk.mat'),'RefWalk','headers');
        load(fullfile(Set.datapath,Set.SubjName,'PertWalk.mat'),'PertWalk');
        Set.headers = headers;
        FileInfo = [];
        TrialInfo = [];


        %% Get all the data
        % unperturbed walking: Get at least n cycles without missing data
        i = 16;
        CombDatSel = [];
        ct = 1;
        while ct < Set.nUnp+1 && i <= length(RefWalk.(Set.WalkSpeed).trial)
            RefSel = RefWalk.(Set.WalkSpeed).trial(i);
            RefSel.iSet = 0;
            [BoolInclude] = ExcludeDataTwente(RefSel,headers,Set.WalkSpeed,s);
            if isempty(RefSel.NaNs) && BoolInclude
                CombDatSel = [CombDatSel RefSel];
                FileInfo = [FileInfo 0];
                TrialInfo = [TrialInfo i];
                ct = ct+1;
            end
            i = i+1;
        end

        % perturbed walking data, two trials per perturbation type, no NaNs in the
        % splines (i.e. missing data).
        for iP = Set.iPertSel
            ct = 1; i= 1;
            while ct<Set.nP+1 && i <= length(PertWalk.(Set.WalkSpeed).PType(iP).trial)
                PSel = PertWalk.(Set.WalkSpeed).PType(iP).trial(i);
                PSel.iSet = iP;
                [BoolInclude] = ExcludeDataTwente(PSel,headers,Set.WalkSpeed,s);
                if isempty(PSel.NaNs) && BoolInclude
                    CombDatSel = [CombDatSel PSel];
                    FileInfo = [FileInfo iP];
                    ct = ct+1;
                    TrialInfo = [TrialInfo i];
                end
                i = i+1;
            end
        end
        Set.FileInfo = FileInfo;

        % get the muscle-tendon information
        [CombDatSel] = GetMuscleModel_vTwente(Set.PolyPath,Set.ModelPath,...
            Set.MuscleNames,CombDatSel,headers,Set.BoolRad,Set.qa_offset);

        % get the Unpeturbed COM trajectories:
        AllRefData = RefWalk.(Set.WalkSpeed);
        RefCycles = GetRefCycle_Twente(AllRefData,Set.headers,Set.BodyMass);

        % Default Geyer model
        Set.ReflexGuess     = [0.02 1.2 0.02 0.4 1 0.2];
        Set.GainsEstimation = 'GeyerDefault';
        Set.MeshFreq        = 1000;
        Set.diaryName = fullfile(Set.ResultsFolder,['Default_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
        [All_GeyerDefault]  = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);

        % Default Geyer model + COM feedback
        Set.ReflexGuess     = [0.02 1.2 0.02 0.4 1 0.2 0.06 -0.1];
        Set.GainsEstimation = 'Geyer_COMd';
        Set.MeshFreq        = 1000;
        Set.diaryName = fullfile(Set.ResultsFolder,['COMfb_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
        [All_GeyerCOM]      = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);

        % Default Geyer model + COMd feedback
        Set.ReflexGuess     = [0.02 1.2 0.02 0.4 1 0.2 0.06 -0.1];
        Set.GainsEstimation = 'Geyer_COMd_GRF';
        Set.MeshFreq        = 1000;
        Set.diaryName       = fullfile(Set.ResultsFolder,['COMfb_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
        [All_GeyerCOM_GRF]  = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);

        % Default Geyer model + COMd feedback
        Set.ReflexGuess     = [1.2 0.02 0.4 1 0.2 0.06 -0.1];
        Set.GainsEstimation = 'Geyer_COMd_GRF_lim';
        Set.MeshFreq        = 1000;
        Set.diaryName       = fullfile(Set.ResultsFolder,['COMfb_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
        Set.Sol.e0          = 0.02;
        [All_GeyerCOM_GRF_lim]  = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);

        %% save results
        save(fullfile(Set.ResultsFolder,[Set.OutName_Gen '_s_' num2str(s) '.mat']),...
            'All_GeyerCOM','All_GeyerDefault','All_GeyerCOM_GRF','FileInfo',...
            'All_GeyerCOM_GRF_lim','TrialInfo');
        clear CombDatSel RefWalk headers PertWalk PID_comb FileInfo
        catch
            % logs all subjects with some errors in the optimization
            Batch.SubjError(s,ib) = 1;
        end
    end
end
save(fullfile(Set.ResultsFolder,'Batch.mat'),'Batch','Set');
