%% Analyse validation on cardiff dataset
%----------------------------------------

%% Analyse prediction medium perturbation

% path information
% MainPath        = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\Geyer_InverseID';
MainPath = 'C:\Users\Maarten\Documents\Software\Sim\Geyer_InverseID';
ResultsFolder   = fullfile(MainPath,'Results\ID_Twente\ID_all_vGRF_TibFix');
figPath         = fullfile(MainPath,'figs','google_docs');

Cols = linspecer(4);
% ColGeyer = [37, 186, 124]./251;
% ColCOMd = [32, 61, 179]./251;
% ColGRF = [176, 18, 52]./251;
ColGeyer = Cols(3,:);
ColGRF = Cols(2,:);

TrackError = nan(7,3,2,2);

Speeds = {'Slow','Fast'};
for sp = 1:2    
    ct = 1;
    for s = [4 5 7 8 10]
        % load the validation simulations (PredictNovelMagnitude_Cardiff.m)
        load(fullfile(ResultsFolder,['Validation_s' num2str(s) '_Speed' Speeds{sp} '.mat']),'ForwardGRF',...
            'ForwardDefault','ForwardGRF_lim');
        % Compute the tracking error
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardDefault,1);
        TrackError(ct,1,1,sp) = Tracking.UnpAll;
        TrackError(ct,1,2,sp) = Tracking.PertAll;
        
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardGRF_lim,1);
        TrackError(ct,2,1,sp) = Tracking.UnpAll;
        TrackError(ct,2,2,sp) = Tracking.PertAll;
        
        Tracking = ComputeTrackingError_ForwardSim_Twente(ForwardGRF,1);
        TrackError(ct,3,1,sp) = Tracking.UnpAll;
        TrackError(ct,3,2,sp) = Tracking.PertAll;
        ct = ct+1;
    end
end

%% plot figure for paper
figure();
subplot(1,2,1)
PlotBar(1,squeeze(TrackError(:,1,2,1)),ColGeyer);
PlotBar(2,squeeze(TrackError(:,3,2,1)),ColGRF);
subplot(1,2,2)
PlotBar(1,squeeze(TrackError(:,1,2,2)),ColGeyer);
PlotBar(2,squeeze(TrackError(:,3,2,2)),ColGRF);

for i=1:2
    subplot(1,2,i)
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'DefaultGeyer','COMd-GRF'});set(gca,'XTick',1:2);
    set(gca,'YTick',0:5:30);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1.5);
    set(gca,'YLim',[0 30]);
end
subplot(1,2,1)
ylabel('RMSE Torque [Nm]');
title('walking 0.62 m/s');
subplot(1,2,2)
% ylabel('RMSE Torque [Nm]');
title('walking 1.25 m/s');

delete_box
set(gcf,'Position',[416         630        1005         206]);
saveas(gcf,fullfile(figPath,'PelvisP_CrossVal_MediumP.png'),'png');
saveas(gcf,fullfile(figPath,'PelvisP_CrossVal_MediumP.svg'),'svg');

