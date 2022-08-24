%% Plot Gains ID all subjects
%---------------------------

clear all; close all; clc;

% plot results generated with script
% ParamEst_TwenteData_TestAllSubj_Shooting.m

%% Path information

MainPath = 'C:\Users\mat950\Documents\Software\Sim\GeyerAnkle_ParamEstimation';

% relative paths
datapath = fullfile(MainPath,'\Results\ResParamID');
DataVlutters = fullfile(MainPath,'Data');

figPath = fullfile(MainPath,'figs','google_docs');
if ~isfolder(figPath)
    mkdir(figPath)
end

%% Layout settings

ColsTemp = linspecer(4);
ColGeyer = ColsTemp(2,:);
ColGRF = ColsTemp(3,:);
ColCOMd = [0.3, 0.3, 0.8];

%% Read all the data
% results from ParamEst_PelvisPerturb_Shooting.m
load(fullfile(datapath,'Batch.mat'));

% Geyer is the default geyer model
% GeyerCOM is the default geyer model with COM veloctity feedback
% GeyerGRF is the geyer model with COM veloctity feedback and smooth
% transition between stance and swing phase reflex gains based on the
% vertical ground reaction force

% pre allocate output matrices
sVect = [4 5 7 8 10];
ReflexGeyer = nan(Batch.N,length(sVect),8);
ReflexCOM = nan(Batch.N,length(sVect),8);
ReflexGRF = nan(Batch.N,length(sVect),8);
ReflexGRF_lim = nan(Batch.N,length(sVect),8);

JGeyer = nan(Batch.N,length(sVect));
JCOM = nan(Batch.N,length(sVect));
JGRF = nan(Batch.N,length(sVect));
JGRF_lim = nan(Batch.N,length(sVect));

TrackError_Geyer = nan(Batch.N,length(sVect));
TrackError_COM = nan(Batch.N,length(sVect));
TrackError_GRF = nan(Batch.N,length(sVect));
TrackError_GRF_lim = nan(Batch.N,length(sVect));

TrackError_Geyer_Unp = nan(Batch.N,length(sVect));
TrackError_COM_Unp = nan(Batch.N,length(sVect));
TrackError_GRF_Unp = nan(Batch.N,length(sVect));
TrackError_GRF_Unp_lim = nan(Batch.N,length(sVect));

% get reflex gains in each batch
for ib = 1:Batch.N
    for s = 1:length(sVect)
        OutName_Gen = Batch.OutNames{ib};
        load(fullfile(datapath,[OutName_Gen '_s_' num2str(sVect(s)) '.mat']),...
            'All_GeyerCOM','All_GeyerDefault','All_GeyerCOM_GRF','All_GeyerCOM_GRF_lim');
        % get the reflex gains
        ReflexGeyer(ib,s,:) =  UnpackAnkleReflex(All_GeyerDefault.Reflex);
        ReflexCOM(ib,s,:) = UnpackAnkleReflex(All_GeyerCOM.Reflex);
        ReflexGRF(ib,s,:) = UnpackAnkleReflex(All_GeyerCOM_GRF.Reflex);
        ReflexGRF_lim(ib,s,:) = UnpackAnkleReflex(All_GeyerCOM_GRF_lim.Reflex);
        % get the cost function
        JGeyer(ib,s) = sum(abs(All_GeyerDefault.J))/length(All_GeyerDefault.J);
        JCOM(ib,s) = sum(abs(All_GeyerCOM.J))/length(All_GeyerCOM.J);
        JGRF(ib,s) = sum(abs(All_GeyerCOM_GRF.J))/length(All_GeyerCOM_GRF.J);
        JGRF_lim(ib,s) = sum(abs(All_GeyerCOM_GRF_lim.J))/length(All_GeyerCOM_GRF_lim.J);
        % get the tracking error
        [TrackError_Geyer_Unp(ib,s), TrackError_Geyer(ib,s), TrackError_Geyer_Push(ib,s), TrackError_Geyer_Pull(ib,s)] = ...
            ComputeTrackingError_Twente(All_GeyerDefault,Batch.nUnpSel{ib});
        [TrackError_COM_Unp(ib,s), TrackError_COM(ib,s), TrackError_COM_Push(ib,s), TrackError_COM_Pull(ib,s)] = ...
            ComputeTrackingError_Twente(All_GeyerCOM,Batch.nUnpSel{ib});
        [TrackError_GRF_Unp(ib,s), TrackError_GRF(ib,s), TrackError_GRF_Push(ib,s), TrackError_GRF_Pull(ib,s)] = ...
            ComputeTrackingError_Twente(All_GeyerCOM_GRF,Batch.nUnpSel{ib});
        [TrackError_GRF_Unp_lim(ib,s), TrackError_GRF_lim(ib,s), TrackError_GRF_lim_Push(ib,s), TrackError_GRF_lim_Pull(ib,s)] = ...
            ComputeTrackingError_Twente(All_GeyerCOM_GRF_lim,Batch.nUnpSel{ib});

    end
end

%% Load validation results
% results from PredictNovelMagnitude.m
TrackError_Val = nan(7,3,2,2);
Speeds = {'Slow'};
for sp = 1:length(Speeds)    
    ct = 1;
    for s = sVect
        % load the validation simulations (PredictNovelMagnitude_Cardiff.m)
        load(fullfile(datapath,['Validation_s' num2str(s) '_Speed' Speeds{sp} '.mat']),'ForwardGRF',...
            'ForwardDefault','ForwardGRF_lim');
        % Compute the tracking error
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardDefault,1);
        TrackError_Val(ct,1,1,sp) = Tracking.UnpAll;
        TrackError_Val(ct,1,2,sp) = Tracking.PertAll;
        
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardGRF_lim,1);
        TrackError_Val(ct,2,1,sp) = Tracking.UnpAll;
        TrackError_Val(ct,2,2,sp) = Tracking.PertAll;
        
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardGRF,1);
        TrackError_Val(ct,3,1,sp) = Tracking.UnpAll;
        TrackError_Val(ct,3,2,sp) = Tracking.PertAll;
        ct = ct+1;
    end
end


%% Display the feedback gains used in the controller

disp('Feedback gains used in controller: ');
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_AllMag');
gains = squeeze(ReflexGRF(iBatch,1,:));
Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};

for i=1:length(Reflex_All_Headers)
    disp(['  ' Reflex_All_Headers{i} '   ' num2str(gains(i))]);
end

%% Figure with tracking error

% plots tracking error in perturbed walking, in perturbed walking and on a
% novel trial not used in the parameter estimation

Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Slow_AllDir_AllMag','Slow_AllDir_AllMag'};
TitlesSel ={'Unperturbed','Perturbed','Validation - Perturbed'};

ColsTemp = linspecer(4);
C1 = ColsTemp(3,:);
C2 = ColsTemp(2,:);

figure('Name','TrackingError Slow Fig Paper');
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;
iSel = iSel_Slow;

subplot(1,3,1)
iBatch = strcmp(Batch.OutNames,NamesSel{1});
PlotBar(1,TrackError_Geyer_Unp(iBatch,iSel),C1);
PlotBar(2,TrackError_GRF_Unp(iBatch,iSel),C2);
set(gca,'XTick',1:2);
% set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});
title('Unperturbed');
set(gca,'FontSize',10);
set(gca,'LineWidth',1.5);
set(gca,'YLim',[0 22]);
set(gca,'YTick',[0 10 20]);

subplot(1,3,2)
iBatch = strcmp(Batch.OutNames,NamesSel{2});
PlotBar(1,TrackError_Geyer(iBatch,iSel),C1);
PlotBar(2,TrackError_GRF(iBatch,iSel),C2);
set(gca,'XTick',1:2);
% set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});
title('Perturbed');
set(gca,'FontSize',10);
set(gca,'LineWidth',1.5);
set(gca,'YLim',[0 22]);
set(gca,'YTick',[0 10 20]);

subplot(1,3,3)
iBatch = strcmp(Batch.OutNames,NamesSel{2});
dsel = squeeze(TrackError_Val(:,1,2,1));
PlotBar(1,dsel,C1);
dsel = squeeze(TrackError_Val(:,3,2,1));
PlotBar(2,dsel,C2);
set(gca,'XTick',1:2);
% set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});
title('Perturbed Validation');
set(gca,'FontSize',10);
set(gca,'LineWidth',1.5);
set(gca,'YLim',[0 22]);
set(gca,'YTick',[0 10 20]);

delete_box
subplot(1,3,1)
ylabel('RMSE Torque [Nm]');
subplot(1,3,2)
ylabel('RMSE Torque [Nm]');
subplot(1,3,3)
ylabel('RMSE Torque [Nm]');

set(gcf,'Position',[896   670   513   202]);
saveas(gcf,fullfile(figPath,'Figure1Paper.png'),'png');
saveas(gcf,fullfile(figPath,'Figure1Paper.svg'),'svg');




%% Summary controller behavior and experiment (joint torque)
%----------------------------------------------------------

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; %
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
%plot figure
HeaderSubPlot = {'Unperturbed','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};
[h1] = PlotTrackingTwente_TorqueAv(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'Figure1_Torque.png'),'png');
saveas(gcf,fullfile(figPath,'Figure1_Torque.svg'),'svg');
saveas(gcf,fullfile(figPath,'Figure1_Torque.fig'),'fig');



%% Appendix: Figure with predicted and measured muscle responses
%----------------------------------------

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
%plot figure
HeaderSubPlot = {'Unperturbed','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};

% load muscle activity - mat file with raw data
Set.SubjName = ['pp_' num2str(SubjSel)];
load(fullfile(DataVlutters,Set.SubjName,'PertWalk.mat'),'PertWalk');
load(fullfile(DataVlutters,Set.SubjName,'RefWalk.mat'),'RefWalk','headers');
% plot activity
[h1] = PlotTrackingTwente_MuscleAv_Exp(datafile,RefWalk,PertWalk,...
    Batch.nUnpSel{iBatch},Batch.nP{iBatch},Batch.iPertSel{iBatch},HeaderSubPlot);


%% Appendix figure to visualise all trials included in parameter estimation
% this figure is not included in the paper, but gives a good overview of
% the different trials that were included in the paper. you can also change
% the subject ID if you want to select you own "representative" subject.
set(0, 'DefaultFigureRenderer', 'opengl');

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; %
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% plot tracking of moment as a functin of time
HeaderSubPlot = {'Unperturbed','Pull 4%','Pull 8%','Pull 12%','Pull 16%',...
   'Push 4%','Push 8%','Push 12%','Push 16%'};
[h1] = PlotTrackingTwente(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);

saveas(gcf,fullfile(figPath,'AllTrialsParamEst.png'),'png');
saveas(gcf,fullfile(figPath,'AllTrialsParamEst.svg'),'svg');

set(0, 'DefaultFigureRenderer', 'painters');


%% Appendix: Smooth transition stance and swing phase

NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 5; %
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% load datafile
D = load(datafile);
% layout
lw = 2;
ColNoGRF = [0.3 0.3 0.8];
% Plot some unperturbed gait cycles
figure();
set(gcf,'Position',[196.3333  708.3333  915.6667  298.6667]);
subplot(1,5,1:3);

% get an unperturbed gait cycle
iUnp = D.FileInfo == 0;
iUnp = iUnp(5);

% plot the inverse dynamic moment
Set = D.All_GeyerCOM.Set;
if isfield(Set,'ScaleTroque') && Set.ScaleTroque
    Tid = D.All_GeyerCOM.Tid./(Set.BodyMass*Set.Height).*(70*1.75);
else
    Tid = D.All_GeyerCOM.Tid;
end
plot(D.All_GeyerCOM.t(iSel(1:end-1)), Tid(iSel(1:end-1)),'--k','LineWidth',lw); hold on;

% plot model simulation with or without smooth transition
iSel = D.All_GeyerCOM.iDatSet(iUnp)+1:D.All_GeyerCOM.iDatSet(iUnp+1);
plot(D.All_GeyerCOM.t(iSel),full(D.All_GeyerCOM.Tmus(iSel)),'Color',ColNoGRF,'LineWidth',lw); hold on;
iSel = D.All_GeyerCOM_GRF.iDatSet(iUnp)+1:D.All_GeyerCOM_GRF.iDatSet(iUnp+1);
plot(D.All_GeyerCOM_GRF.t(iSel),full(D.All_GeyerCOM_GRF.Tmus(iSel)),'Color',ColGeyer,'LineWidth',lw); hold on;
legend('Inverse dynamics', 'Model discrete transition gains', 'Model smooth transition');
xlabel('Time [s]');
ylabel('Ankle moment [Nm]');

% plot the average tracking error in all subjects
subplot(1,5,4)
PlotBar(4,TrackError_COM_Unp(iBatch,:),ColNoGRF); hold on;
PlotBar(5, TrackError_GRF_Unp_lim(iBatch,:),ColGeyer);
ylabel('RMSE [Nm]')
title('Unperturbed')
subplot(1,5,5)
PlotBar(1,TrackError_COM(iBatch,:),ColNoGRF); hold on;
PlotBar(2, TrackError_GRF_lim(iBatch,:),ColGeyer)
ylabel('RMSE [Nm]')
title('Perturbed')

% adapt the 
iPlt = {1:3,4,5};
for i=1:length(iPlt)
    subplot(1,5,iPlt{i})
    set(gca,'box','off');
    set(gca,'LineWidth',1.4);
    set(gca,'FontSize',12);
end

%% Appendix Added model without smooth transition in the overview figure
% this figure is not included in the paper, but gives a good overview of
% the different trials that were included in the paper. you can also change
% the subject ID if you want to select you own "representative" subject.
set(0, 'DefaultFigureRenderer', 'opengl');

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 5; %
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% plot tracking of moment as a functin of time
HeaderSubPlot = {'Unperturbed','Pull 4%','Pull 8%','Pull 12%','Pull 16%',...
   'Push 4%','Push 8%','Push 12%','Push 16%'};
[h1] = PlotTrackingTwente_vSmoothTransition(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);

%% Appendix figure direction independent gains

NamesSel = {'Slow_AllDir_AllMag','Slow_AllDir_Mag34','Slow_AllDir_Mag78'};
figure();

subplot(1,2,1)
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_AllMag');
PlotBar(1,TrackError_GRF_Pull(iBatch,:),ColGRF);
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_Mag34');
PlotBar(2,TrackError_GRF_Pull(iBatch,:),[0.6 0.6 0.6]);
title('pull perturbations')

subplot(1,2,2)
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_AllMag');
PlotBar(1,TrackError_GRF_Push(iBatch,:),ColGRF);
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_Mag78');
PlotBar(2,TrackError_GRF_Push(iBatch,:),[0.6 0.6 0.6]);
title('push perturbations')

for i=1:2
    subplot(1,2,i)
    ylabel('RMSE []');
    set(gca,'box','off');
    set(gca,'LineWidth',1.4);
    set(gca,'FontSize',12);
    set(gca,'XTick',[]);
end
set(gcf,'Position',[275.0000  802.3333  837.0000  204.6667]);
saveas(gcf,fullfile(figPath,'AppendixDirDepGains.svg'),'svg');
