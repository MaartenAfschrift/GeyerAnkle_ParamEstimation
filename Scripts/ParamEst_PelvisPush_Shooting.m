
%% Gain estimation data Twente
%-----------------------------------

clear all; close all; clc;

% path information
Set.datapath = 'C:\Users\mat950\Documents\Data\Vlutters2018\ParamID';
Set.MainPath = 'C:\Users\mat950\Documents\Software\Sim\GeyerAnkle_ParamEstimation';
Set.PolyPath        = fullfile(Set.MainPath,'OsimModel\MuscleAnalysis');
% Set.MuscleInfoPath  = fullfile(Set.MainPath,'Functions\HillModel');
Set.ModelPath       = fullfile(Set.MainPath,'OsimModel\Models','gait_2D_Strength.osim');
Set.ResultsFolder   = fullfile(Set.MainPath,'Results','Example');

% add paths
addpath(genpath(fullfile(Set.MainPath,'Functions')));
import casadi.*

% Select Dataset
Set.TwenteData = 1;   % Boolean for selecting data Twente

% subject IDs5
Set.SubjIDs = 5; %               Set.SubjIDs = [1 4 5 6 7 8 10];

% Body mass
Set.BodyMass = 70;          % in kg

% angles in radians or deg to evaluate polynomials
Set.BoolRad = false;

% set info from the model
Set.MuscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};

% initial guess based on previous solution
Set.PathIG = fullfile(pwd,'Reflex_AddHip.mat');

% boolean to plot figures
Set.BoolPlot = 0;

% scaling of the feedback gains
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

% Vlugt
Set.Vlugt2010 = false;
Set.PlotVlugt = false;

% constant "band around COM data"
Set.ConstantSTD = true;
Set.COMref_std = 0.04; % 2 cm deviation from average
Set.COMdref_std = 0.06; % 6 cm/s deviation in COM velocity from average
Set.COMddref_std = 0.01; % 1cm/s2 deviation in COM acceleration from average

% perturbation settings
Set.COMfb = 1;
Set.COM_std = false;
Set.PercGaitCycle = true;
Set.COM_TimeDelay = 0.03;

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


%% Create a Batch for simulations

Batch.Speeds ={'Slow','Slow','Slow','Fast','Fast','Fast',...
    'Slow','Slow','Slow','Slow',...
    'Fast','Fast','Fast','Fast',...
    'Slow','Fast'};
% 1->4: backward increasing magnitude   5->8: forward increasing magnitude
Batch.iPertSel={1:8,1:4,5:8,1:8,1:4,5:8,...
    1:2,3:4,6:7,7:8,...
    1:2,3:4,6:7,7:8,...
    [1 2 4 5 6 8],[1 2 4 5 6 8]}; 
Batch.nUnpSel = {8,8,8,8,8,8,...
    4,4,4,4,...
    4,4,4,4,...
    8,8}; % number of unperturbed gait cycles
Batch.nP = {2,2,2,2,2,2,...
    3,3,3,3,...
    3,3,3,3,...
    2,2};      % number of perturbed gait cycles
Batch.N = length(Batch.Speeds);                   % number of batches
Batch.OutNames = {'Slow_AllDir_AllMag','Slow_Backward_AllMag','Slow_Forward_AllMag',...
    'Fast_AllDir_AllMag','Fast_Backward_AllMag','Fast_Forward_AllMag',...
    'Slow_AllDir_Mag12','Slow_AllDir_Mag34','Slow_AllDir_Mag56','Slow_AllDir_Mag78',...
    'Fast_AllDir_Mag12','Fast_AllDir_Mag34','Fast_AllDir_Mag56','Fast_AllDir_Mag78',...
    'Slow_Valid_Excl37','Fast_Valid_Excl37'};
Batch.SubjError = zeros(10,Batch.N);
SetO = Set;

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
%         try
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
            
            % Default Geyer model + COMd and COMdd feedback
            Set.ReflexGuess     = [0.02 1.4 0.02 0.4 1 0.2 0 0 0.3 -0.2 0.1 0];
            Set.GainsEstimation = 'Geyer_COMAll_GRF';
            Set.MeshFreq        = 1000;
            Set.diaryName       = fullfile(Set.ResultsFolder,['COMfb_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
            [All_COMAll]  = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);
            
            % Default Geyer model + COMd and COMdd feedback
            Set.ReflexGuess     = [0.02 1.4 0.02 0.4 1 0.2 0.3 -0.2 0.1 0];
            Set.GainsEstimation = 'Geyer_COMdd_GRF';
            Set.MeshFreq        = 1000;
            Set.diaryName       = fullfile(Set.ResultsFolder,['COMfb_' Set.OutName_Gen '_s_' num2str(s) '.txt']);
            [All_COMdd]  = GainEstimation_AnkleOnly_nCycl_Shooting_vTwente(CombDatSel,Set,RefCycles);
            
            %% Plot Main Results
%             h1 = figure();
% %             set(h1,'Position',get(0,'ScreenSize'));
%             set(h1,'Position',[680          89         788        1191]);
%             for i=1:length(CombDatSel)
%                 subplot(ceil(length(CombDatSel)/2),2,i);
%                 
%                 % default geyer model shooting
%                 iSel = All_GeyerDefault.iDatSet(i)+1:All_GeyerDefault.iDatSet(i+1);
%                 plot(All_GeyerDefault.t(iSel),full(All_GeyerDefault.Tmus(iSel)),'r'); hold on;
%                 
%                 % COM geyer model shooting
%                 iSel = All_GeyerCOM.iDatSet(i)+1:All_GeyerCOM.iDatSet(i+1);
%                 plot(All_GeyerCOM.t(iSel),full(All_GeyerCOM.Tmus(iSel)),'b'); hold on;
%                 
%                 % COM geyer GRF model shooting
%                 iSel = All_GeyerCOM_GRF.iDatSet(i)+1:All_GeyerCOM_GRF.iDatSet(i+1);
%                 plot(All_GeyerCOM_GRF.t(iSel),full(All_GeyerCOM_GRF.Tmus(iSel)),'g'); hold on;
%                 
%                 % ID moment
%                 plot(All_GeyerCOM.t(iSel), All_GeyerCOM.Tid(iSel),'--k');
%                 
%                 ylabel('Ankle moment [Nm]');
%                 set(gca,'YLim',[-100 30]);
%             end
%             legend('DefaultGeyer','Geyer COM','Geyer COM GRF','Tid');
% 
%             saveas(h1,fullfile(Set.ResultsFolder,[Set.OutName_Gen '_s_' num2str(s) '.fig']));
%             saveas(h1,fullfile(Set.ResultsFolder,[Set.OutName_Gen '_s_' num2str(s) '.svg']),'svg');
%             saveas(h1,fullfile(Set.ResultsFolder,[Set.OutName_Gen '_s_' num2str(s) '.png']),'png');
            
            %% save results
            save(fullfile(Set.ResultsFolder,[Set.OutName_Gen '_s_' num2str(s) '.mat']),...
                'All_GeyerCOM','All_GeyerDefault','All_GeyerCOM_GRF','FileInfo',...
                'All_GeyerCOM_GRF_lim','TrialInfo','All_COMAll','All_COMdd');
            clear CombDatSel RefWalk headers PertWalk PID_comb FileInfo
%         catch
%             Batch.SubjError(s,ib) = 1;
%         end
    end
end
save(fullfile(Set.ResultsFolder,'Batch.mat'),'Batch','Set');
