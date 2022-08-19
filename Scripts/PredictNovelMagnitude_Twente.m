%% Predict novel magnitude Twente data
%-------------------------------------

% this script is used to run a forward simulation on a perturbation trials
% that was not included in the paramter estimation process

clear all; clc;
% path information
MainPath        = 'C:\Users\Maarten\Documents\Software\Sim\GeyerAnkle_ParamEstimation';
resultfolder    = fullfile(MainPath,'Results','ID_Twente','ID_COMdd_Fix_v2');
datapath        = 'D:\DataGeyerID';    % path with data of Vlutters 2018 organized in format for parameter estimation


PolyPath        = fullfile(MainPath,'OsimModel\MuscleAnalysis');
ModelPath       = fullfile(MainPath,'OsimModel\Models','gait_2D_Strength.osim');


% load the batch (structure with all simulations)
load(fullfile(resultfolder,'Batch.mat'),'Batch','Set');

Validations = {'Slow_Valid_Excl37','Fast_Valid_Excl37'};


% Run forward simulation
%% Collect data to predict perturbation magnitude 3 and 7
for ibn = 1:length(Validations)
    ib = strcmp(Batch.OutNames,Validations{ibn});
    Set.OutName_Gen = Batch.OutNames{ib};
    for s = [4 5 7 8 10]
        
        % load the optimized parameters
        load(fullfile(resultfolder,[Set.OutName_Gen '_s_' num2str(s) '.mat']),...
                'All_GeyerCOM','All_GeyerDefault','All_GeyerCOM_GRF','FileInfo',...
                'All_GeyerCOM_GRF_lim','TrialInfo');
        
        Set = All_GeyerCOM_GRF_lim.Set;
        Set.OutName_Gen = Batch.OutNames{ib};
        Set.WalkSpeed = Batch.Speeds{ib}; % walking speed
        Set.iPertSel = [3 7];
        Set.nUnp = 0;
        Set.nP = 3;
           
        
        %% Load the data of this subject
        Set.SubjName = ['pp_' num2str(Set.sID)];
        load(fullfile(datapath,Set.SubjName,'RefWalk.mat'),'RefWalk','headers');
        load(fullfile(datapath,Set.SubjName,'PertWalk.mat'),'PertWalk');
        Set.headers = headers;
        
        %% get refcycles
        AllRefData = RefWalk.(Set.WalkSpeed);
        RefCycles = GetRefCycle_Twente(AllRefData,Set.headers);
       
        
        %% Get all the data
        % unperturbed walking: Get at least n cycles without missing data
        i = max(TrialInfo)+1;
        CombDatSel = [];
        ct = 1;
        FileInfo2 = [];
        TrialInfo2 = [];
        while ct < Set.nUnp+1 && i <= length(RefWalk.(Set.WalkSpeed).trial)
            RefSel = RefWalk.(Set.WalkSpeed).trial(i);
            [BoolInclude] = ExcludeDataTwente(RefSel,headers,Set.WalkSpeed,s);
            if isempty(RefSel.NaNs) && BoolInclude
                CombDatSel = [CombDatSel RefSel];
                FileInfo2 = [FileInfo2 0];
                TrialInfo2 = [TrialInfo2 i];
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
                [BoolInclude] = ExcludeDataTwente(PSel,headers,Set.WalkSpeed,s);
                if isempty(PSel.NaNs) && BoolInclude
                    CombDatSel = [CombDatSel PSel];
                    FileInfo2 = [FileInfo2 iP];
                    ct = ct+1;
                    TrialInfo2 = [TrialInfo2 i];
                end
                i = i+1;
            end
        end
        
        % get the muscle-tendon information
        [CombDatSel] = GetMuscleModel_vTwente(PolyPath,ModelPath,...
            Set.MuscleNames,CombDatSel,headers,Set.BoolRad,Set.qa_offset);
        
        % Run forward simulation on selected data
        for i=1:length(CombDatSel)
            ForwardGRF(i) = ForwardSim_GeyerModel_Explicit_AnkleOnly_vTwente(CombDatSel(i),...
                All_GeyerCOM_GRF.Set,All_GeyerCOM_GRF.Reflex,RefCycles);
        end
        
        for i=1:length(CombDatSel)
            ForwardGRF_lim(i) = ForwardSim_GeyerModel_Explicit_AnkleOnly_vTwente(CombDatSel(i),...
                All_GeyerCOM_GRF_lim.Set,All_GeyerCOM_GRF_lim.Reflex,RefCycles);
        end
        
        for i=1:length(CombDatSel)
            ForwardDefault(i) = ForwardSim_GeyerModel_Explicit_AnkleOnly_vTwente(CombDatSel(i),...
                All_GeyerDefault.Set,All_GeyerDefault.Reflex,RefCycles);
        end
        
        figure();
        for i=1:length(ForwardGRF)
            subplot(length(ForwardGRF),1,i)
            plot(ForwardGRF(i).t,ForwardGRF(i).Tankle,'r'); hold on;
            plot(ForwardDefault(i).t,ForwardDefault(i).Tankle,'b'); hold on;
            plot(ForwardDefault(i).t,ForwardDefault(i).Tid_ankle,'--k'); hold on;
        end
%         disp('test');
        save(fullfile(resultfolder,['Validation_s' num2str(s) '_Speed' Set.WalkSpeed '.mat']),'ForwardGRF',...
            'ForwardDefault','ForwardGRF_lim');
        clear ForwardGRF ForwardDefault ForwardGRF_lim
    end
end
