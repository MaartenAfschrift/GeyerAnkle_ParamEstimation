%% Plot Gains ID all subjects
%---------------------------

clear all; close all; clc;

% plot results generated with script
% ParamEst_TwenteData_TestAllSubj_Shooting.m

%% Path information

MainPath = 'C:\Users\Maarten\Documents\Software\Sim\GeyerAnkle_ParamEstimation';
datapath = fullfile(MainPath,'\Results\ID_Twente\ID_COMdd_Fix_v2');
DataVlutters = 'D:\DataVlutters';

figPath = fullfile(MainPath,'figs','google_docs');
if ~isfolder(figPath)
    mkdir(figPath)
end

%% Layout settings
ColGeyer = [37, 186, 124]./251;
ColCOMd = [32, 61, 179]./251;
ColGRF = [176, 18, 52]./251;

%% Read all the data
load(fullfile(datapath,'Batch.mat'));

% sVect = [1 4 5 6 7 8 10];
sVect = [4 5 7 8 10];
ReflexGeyer = nan(Batch.N,length(sVect),8);
ReflexCOM = nan(Batch.N,length(sVect),8);
ReflexGRF = nan(Batch.N,length(sVect),8);
JGeyer = nan(Batch.N,length(sVect));
JCOM = nan(Batch.N,length(sVect));
JGRF = nan(Batch.N,length(sVect));
TrackError_Geyer = nan(Batch.N,length(sVect));
TrackError_COM = nan(Batch.N,length(sVect));
TrackError_GRF = nan(Batch.N,length(sVect));
TrackError_Geyer_Unp = nan(Batch.N,length(sVect));
TrackError_COM_Unp = nan(Batch.N,length(sVect));
TrackError_GRF_Unp = nan(Batch.N,length(sVect));


% get reflex gains in each batch
for ib = 1:Batch.N
    for s = 1:length(sVect)
%         try
            OutName_Gen = Batch.OutNames{ib};
            load(fullfile(datapath,[OutName_Gen '_s_' num2str(sVect(s)) '.mat']),...
                'All_GeyerCOM','All_GeyerDefault','All_GeyerCOM_GRF_lim');
            ReflexGeyer(ib,s,:) =  UnpackAnkleReflex(All_GeyerDefault.Reflex);
            ReflexCOM(ib,s,:) = UnpackAnkleReflex(All_GeyerCOM.Reflex);
            ReflexGRF(ib,s,:) = UnpackAnkleReflex(All_GeyerCOM_GRF_lim.Reflex);
            JGeyer(ib,s) = sum(abs(All_GeyerDefault.J))/length(All_GeyerDefault.J);
            JCOM(ib,s) = sum(abs(All_GeyerCOM.J))/length(All_GeyerCOM.J);
            JGRF(ib,s) = sum(abs(All_GeyerCOM_GRF_lim.J))/length(All_GeyerCOM_GRF_lim.J);
            [TrackError_Geyer_Unp(ib,s), TrackError_Geyer(ib,s)] = ComputeTrackingError_Twente(All_GeyerDefault,Batch.nUnpSel{ib});
            [TrackError_COM_Unp(ib,s), TrackError_COM(ib,s)] = ComputeTrackingError_Twente(All_GeyerCOM,Batch.nUnpSel{ib});
            [TrackError_GRF_Unp(ib,s), TrackError_GRF(ib,s)] = ComputeTrackingError_Twente(All_GeyerCOM_GRF_lim,Batch.nUnpSel{ib});
%         catch
%         end
    end
end

%% Load validation results
TrackError_Val = nan(7,3,2,2);

Speeds = {'Slow','Fast'};
for sp = 1:2    
    ct = 1;
    for s = [4 5 7 8 10]
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

%% Figure with feedback gains


% Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
% Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];
% 
% NamesSel = {'Slow_Backward_AllMag','Fast_Backward_AllMag','Slow_Forward_AllMag','Fast_Forward_AllMag'};
% TitlesSel ={'Slow backward','Fast backward','Slow forward','Fast forward'};
% 
% figure();
% for j = 1:4
%     subplot(2,2,j)
%     iBatch = strcmp(Batch.OutNames,NamesSel{j});
%     for i=1:8
%         PlotBar(i*4-3,squeeze(ReflexGeyer(iBatch,:,i)),ColGeyer);
%         PlotBar(i*4-2,squeeze(ReflexCOM(iBatch,:,i)),ColCOMd);
%         PlotBar(i*4-1,squeeze(ReflexGRF(iBatch,:,i)),ColGRF);
%     end
%     set(gca,'XTick',(1:8)*4-2);
%     set(gca,'XTickLabel',Reflex_All_Headers);
%     set(gca,'XTickLabelRotation',60);
%     title(TitlesSel{j});
%     set(gca,'FontSize',10);
%     set(gca,'LineWidth',1.5);
% end
% delete_box
% saveas(gcf,fullfile(figPath,'PelivPush_Gains.png'),'png');


%% Figure with feedback gains


Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Fast_AllDir_AllMag'};
TitlesSel ={'0.62 m/s','1.25 m/s'};
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;

figure();
for j = 1:2
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    if j ==1
        iSel = iSel_Slow;
    elseif j == 2
        iSel = iSel_Fast;
    end
    for i=1:8
        subplot(2,4,i);
        PlotBar(j,squeeze(ReflexGRF(iBatch,iSel,i)),ColGRF); hold on;
    end
    delete_box
    set(gcf,'Position',[610         584        1066         428]);
end

for i=1:8
    subplot(2,4,i)
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    title(Reflex_All_Headers{i});
    set(gca,'XTick',[1:2]);
    set(gca,'XTickLabel',TitlesSel);
    set(gca,'XTickLabelRotation',60);
end
saveas(gcf,fullfile(figPath,'PelivPush_Gains.png'),'png');

%% Display the feedback gains used in the controller

disp('Feedback gains used in controller: ');
iBatch = strcmp(Batch.OutNames,'Slow_AllDir_AllMag');
gains = squeeze(ReflexGRF(iBatch,2,:));
Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};

for i=1:length(Reflex_All_Headers)
    disp(['  ' Reflex_All_Headers{i} '   ' num2str(gains(i))]);
end
% 
% 
%% Figure with Tracking objective
% 
% 
% Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
% Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];
% 
% NamesSel = {'Slow_Backward_AllMag','Fast_Backward_AllMag','Slow_Forward_AllMag','Fast_Forward_AllMag'};
% TitlesSel ={'Slow backward','Fast backward','Slow forward','Fast forward'};
% 
% figure();
% for j = 1:4
%     subplot(2,2,j)
%     iBatch = strcmp(Batch.OutNames,NamesSel{j});
%     PlotBar(1,JGeyer(iBatch,:)*50,[1 0.2 0.2]);
%     PlotBar(2,JCOM(iBatch,:)*50,[0.2 0.2 1]);
%     set(gca,'XTick',1:2);
%     set(gca,'XTickLabel',{'Default','COM-fb'});
%     set(gca,'XTickLabelRotation',60);
%     title(TitlesSel{j});
%     set(gca,'FontSize',10);
%     set(gca,'LineWidth',1.5);
% end
% delete_box


%% Figure with Tracking error

Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];
NamesSel = {'Slow_Backward_AllMag','Fast_Backward_AllMag','Slow_Forward_AllMag','Fast_Forward_AllMag'};
TitlesSel ={'Slow backward','Fast backward','Slow forward','Fast forward'};
NamesSelAllDir = {'Slow_AllDir_AllMag','Fast_AllDir_AllMag','Slow_AllDir_AllMag','Fast_AllDir_AllMag'};

figure('Name','TrackingError All mag Dir Dep');
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;
for j = 1:4
    subplot(2,2,j)
    if j ==1 || j==3
        iSel = iSel_Slow;
    elseif j == 2 || j==4
        iSel = iSel_Fast;
    end
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer(iBatch,iSel),ColGeyer);
    PlotBar(2,TrackError_COM(iBatch,iSel),ColCOMd);
    PlotBar(3,TrackError_GRF(iBatch,iSel),ColGRF);
    iBatch2 = strcmp(Batch.OutNames,NamesSelAllDir{j});
    PlotBar(4,TrackError_GRF(iBatch2,iSel),[0.6 0.6 0.6]);
    set(gca,'XTick',1:4);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF','GRF-AllDir'});
    set(gca,'XTickLabelRotation',60);
    title(TitlesSel{j});
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 20]);
end
delete_box
for i=1:4
    subplot(2,2,i)
    ylabel('RMSE torque [Nm]');
end
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_DirDependent.png'),'png');


%% Figure with Tracking error


Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Fast_AllDir_AllMag'};
TitlesSel ={'walk 0.62 m/s','walk 1.25 m/s'};

figure('Name','TrackingError All P');
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;
for j = 1:2
    if j ==1
        iSel = iSel_Slow;
    elseif j == 2
        iSel = iSel_Fast;
    end
    subplot(1,2,j)
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer(iBatch,iSel),ColGeyer);
    PlotBar(2,TrackError_COM(iBatch,iSel),ColCOMd);
    PlotBar(3,TrackError_GRF(iBatch,iSel),ColGRF);
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF'});
    title(TitlesSel{j});
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 22]);
end
delete_box
subplot(1,2,1)
ylabel('RMSE Torque [Nm]');
% leg1 = legend({'DefaultGeyer','COMd','COMd-GRF'});
%
% set(leg1,...
%     'Position',[0.367272421755037 0.964852672238066 0.285100282391368 0.0249632887115563],...
%     'Orientation','horizontal');
set(gcf,'Position',[916   808   969   299]);
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag.png'),'png');
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag.svg'),'svg');

%% Figure with Tracking error with Unp as well



Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Fast_AllDir_AllMag'};
TitlesSel ={'0.62 m/s','1.25 m/s'};

figure('Name','TrackingError All P');
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;
for j = 1:2
    if j ==1
        iSel = iSel_Slow;
    elseif j == 2
        iSel = iSel_Fast;
    end
    subplot(2,2,j)
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer_Unp(iBatch,iSel),ColGeyer);
    PlotBar(2,TrackError_COM_Unp(iBatch,iSel),ColCOMd);
    PlotBar(3,TrackError_GRF_Unp(iBatch,iSel),ColGRF);
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF'});
    title(['Unperturbed ' TitlesSel{j}]);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 22]);
    subplot(2,2,j+2)
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer(iBatch,iSel),ColGeyer);
    PlotBar(2,TrackError_COM(iBatch,iSel),ColCOMd);
    PlotBar(3,TrackError_GRF(iBatch,iSel),ColGRF);
    set(gca,'XTick',1:3);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF'});
    title(['Perturbed ' TitlesSel{j}]);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 22]);
end
delete_box
subplot(2,2,1)
ylabel('RMSE Torque [Nm]');
subplot(2,2,3)
ylabel('RMSE Torque [Nm]');
% leg1 = legend({'DefaultGeyer','COMd','COMd-GRF'});
%
% set(leg1,...
%     'Position',[0.367272421755037 0.964852672238066 0.285100282391368 0.0249632887115563],...
%     'Orientation','horizontal');
set(gcf,'Position',[401         667        1290         446]);
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag_UnpPert.png'),'png');

%% Figure with Tracking error with Unp as well



Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Fast_AllDir_AllMag'};
TitlesSel ={'0.62 m/s','1.25 m/s'};

ColsTemp = linspecer(4);
C1 = ColsTemp(3,:);
C2 = ColsTemp(2,:);

figure('Name','TrackingError All P');
iSel_Fast = 1:5;%[2 3 5 6 7]; % Convergence errors in s1 and s6
iSel_Slow = 1:5;
for j = 1:2
    if j ==1
        iSel = iSel_Slow;
    elseif j == 2
        iSel = iSel_Fast;
    end
    subplot(2,2,j)
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer_Unp(iBatch,iSel),C1);
    PlotBar(2,TrackError_GRF_Unp(iBatch,iSel),C2);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});
    title(['Unperturbed ' TitlesSel{j}]);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 22]);
    subplot(2,2,j+2)
    iBatch = strcmp(Batch.OutNames,NamesSel{j});
    PlotBar(1,TrackError_Geyer(iBatch,iSel),C1);
    PlotBar(2,TrackError_GRF(iBatch,iSel),C2);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});
    title(['Perturbed ' TitlesSel{j}]);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 22]);
end
delete_box
subplot(2,2,1)
ylabel('RMSE Torque [Nm]');
subplot(2,2,3)
ylabel('RMSE Torque [Nm]');
% leg1 = legend({'DefaultGeyer','COMd','COMd-GRF'});
%
% set(leg1,...
%     'Position',[0.367272421755037 0.964852672238066 0.285100282391368 0.0249632887115563],...
%     'Orientation','horizontal');
set(gcf,'Position',[492   693   398   425]);
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag_UnpPert_v2.png'),'png');
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag_UnpPert_v2.svg'),'svg');

%% Figure RMSE slow walking paper


Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];

NamesSel = {'Slow_AllDir_AllMag','Slow_AllDir_AllMag','Slow_AllDir_AllMag'};
TitlesSel ={'Unperturbed','Perturbed','Validation'};

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
% leg1 = legend({'DefaultGeyer','COMd','COMd-GRF'});

subplot(1,3,3)
ylabel('RMSE Torque [Nm]');
% leg1 = legend({'DefaultGeyer','COMd','COMd-GRF'});

% set(leg1,...
%     'Position',[0.367272421755037 0.964852672238066 0.285100282391368 0.0249632887115563],...
%     'Orientation','horizontal');
set(gcf,'Position',[896   670   513   202]);
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag_SlowPaper.png'),'png');
saveas(gcf,fullfile(figPath,'PelivPush_TrackingError_AllMag_SlowPaper.svg'),'svg');




%% Figure with Tracking error


% Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
% Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];
% 
% NamesSel = {'Slow_AllDir_Mag12','Fast_AllDir_Mag12','Slow_AllDir_Mag56','Fast_AllDir_Mag56'};
% TitlesSel ={'Slow backward SmallP','Fast backward SmallP','Slow forward SmallP','Fast forward SmallP'};
% 
% figure('Name','TrackingError Small perturbations');
% for j = 1:4
%     subplot(2,2,j)
%     iBatch = strcmp(Batch.OutNames,NamesSel{j});
%     PlotBar(1,TrackError_Geyer(iBatch,:),ColGeyer);
%     PlotBar(2,TrackError_COM(iBatch,:),ColCOMd);
%     PlotBar(3,TrackError_GRF(iBatch,:),ColGRF);
%     set(gca,'XTick',1:3);
%     set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF'});
%     title(TitlesSel{j});
%     set(gca,'FontSize',10);
%     set(gca,'LineWidth',1.5);
% end
% delete_box

%% Figure with Tracking error
% Reflex_All_Headers = {'Sol e0','Sol G','Tib e0', 'Tib G','Tib loff','Tib Gsol'};
% Reflex_All_Headers = [Reflex_All_Headers {'Sol COMd','Tib COMd'}];
% 
% NamesSel = {'Slow_AllDir_Mag34','Fast_AllDir_Mag34','Slow_AllDir_Mag78','Fast_AllDir_Mag78'};
% TitlesSel ={'Slow backward LargeP','Fast backward SmallP','Slow forward LargeP','Fast forward LargeP'};
% 
% figure('Name','TrackingError Large perturbations');
% for j = 1:4
%     subplot(2,2,j)
%     iBatch = strcmp(Batch.OutNames,NamesSel{j});
%     PlotBar(1,TrackError_Geyer(iBatch,:),ColGeyer);
%     PlotBar(2,TrackError_COM(iBatch,:),ColCOMd);
%     PlotBar(3,TrackError_GRF(iBatch,:),ColGRF);
%     set(gca,'XTick',1:3);
%     set(gca,'XTickLabel',{'DefaultGeyer','COMd','COMd-GRF'});
%     set(gca,'XTickLabelRotation',60);
%     title(TitlesSel{j});
%     set(gca,'FontSize',10);
%     set(gca,'LineWidth',1.5);
% end
% delete_box

%% Figure with tracking of joint moment as a function of time

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% plot tracking of moment as a functin of time
HeaderSubPlot = {'Unperturbed walking at 0.62 m/s','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};
[h1] = PlotTrackingTwente(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);

saveas(gcf,fullfile(figPath,'PelivPush_TorqueVSTime_Slow.png'),'png');
saveas(gcf,fullfile(figPath,'PelivPush_TorqueVSTime_Slow.svg'),'svg');


%% Figure with tracking of joint moment as a function of time

% name of optimization
NameBatchSel = 'Fast_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% plot tracking of moment as a functin of time
HeaderSubPlot = {'Unperturbed walking at 1.25 m/s','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};
[h1] = PlotTrackingTwente(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'PelivPush_TorqueVSTime_Fast.png'),'png');

%% Figure with muscle resposne

% name of optimization
NameBatchSel = 'Slow_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
% plot tracking of moment as a functin of time
HeaderSubPlot = {'Unperturbed walking at 0.62 m/s','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};
[h1] = PlotTrackingTwente_Soleus(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);

[h2] = PlotTrackingTwente_Tibialis(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
% saveas(gcf,fullfile(figPath,'PelivPush_TorqueVSTime_Fast.png'),'png');


%% Less details in figure muscle response
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

% [h1] = PlotTrackingTwente_MuscleAv(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
%     Batch.iPertSel{iBatch},HeaderSubPlot);
% saveas(gcf,fullfile(figPath,'PelivPush_Slow_MuscleResponse.png'),'png');

% load muscle activity - mat file with raw data
Set.datapath = [DataVlutters '\DataFiles_HSAdapt'];
Set.SubjName = ['pp_' num2str(SubjSel)];
load(fullfile(Set.datapath,Set.SubjName,'PertWalk.mat'),'PertWalk');
load(fullfile(Set.datapath,Set.SubjName,'RefWalk.mat'),'RefWalk','headers');
% plot activity
[h1] = PlotTrackingTwente_MuscleAv_Exp(datafile,RefWalk,PertWalk,...
    Batch.nUnpSel{iBatch},Batch.nP{iBatch},Batch.iPertSel{iBatch},HeaderSubPlot);
% saveas(gcf,fullfile(figPath,'PelivPush_Slow_MuscleResponse_EMG.png'),'png');

%% Less details in figure muscle response
%----------------------------------------

% name of optimization
NameBatchSel = 'Fast_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
%plot figure
HeaderSubPlot = {'Unperturbed','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};

% [h1] = PlotTrackingTwente_MuscleAv(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
%     Batch.iPertSel{iBatch},HeaderSubPlot);
% saveas(gcf,fullfile(figPath,'PelivPush_Fast_MuscleResponse.png'),'png');

% load muscle activity - mat file with raw data
Set.datapath = [DataVlutters '\DataFiles_HSAdapt'];
Set.SubjName = ['pp_' num2str(SubjSel)];
load(fullfile(Set.datapath,Set.SubjName,'PertWalk.mat'),'PertWalk');
load(fullfile(Set.datapath,Set.SubjName,'RefWalk.mat'),'RefWalk','headers');
% plot activity
[h1] = PlotTrackingTwente_MuscleAv_Exp(datafile,RefWalk,PertWalk,...
    Batch.nUnpSel{iBatch},Batch.nP{iBatch},Batch.iPertSel{iBatch},HeaderSubPlot);
% saveas(gcf,fullfile(figPath,'PelivPush_Fast_MuscleResponse_EMG.png'),'png');

%% Less detailed figure in torque response
%-----------------------------------------

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
[h1] = PlotTrackingTwente_TorqueAv(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'PelivPush_Slow_TorqueResponseModels.png'),'png');
saveas(gcf,fullfile(figPath,'PelivPush_Slow_TorqueResponseModels.svg'),'svg');
saveas(gcf,fullfile(figPath,'PelivPush_Slow_TorqueResponseModels.fig'),'fig');


% name of optimization
NameBatchSel = 'Fast_AllDir_AllMag';
% select index in batch structure
iBatch = strcmp(Batch.OutNames,NameBatchSel); 
% select subject
SubjSel = 4; % [1 4 5 6 7 8 10];
% path to selected datafile
datafile = fullfile(datapath,[NameBatchSel '_s_' num2str(SubjSel) '.mat']);
%plot figure
HeaderSubPlot = {'Unperturbed','Pull 1','Pull 2','Pull 3','Pull 4',...
   'Push 1','Push 2','Push 3','Push 4'};
[h1] = PlotTrackingTwente_TorqueAv(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'PelivPush_Fast_TorqueResponseModels.png'),'png');

%% Plot exoskeleton assistance in option IIb-3

% Provide difference between default Geyer and COM geyer as balance support
% assistance.

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
[h1] = PlotTrackingTwente_TorqueAv_Exo(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'PelivPush_Slow_ExoController.png'),'png');

%% Plot exoskeleton assistance Full case

% Provide difference between default Geyer and COM geyer as balance support
% assistance.

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
[h1] = PlotTrackingTwente_TorqueAv_ExoFull(datafile,Batch.nUnpSel{iBatch},Batch.nP{iBatch},...
    Batch.iPertSel{iBatch},HeaderSubPlot);
saveas(gcf,fullfile(figPath,'PelivPush_Slow_ExoControllerFull.png'),'png');


