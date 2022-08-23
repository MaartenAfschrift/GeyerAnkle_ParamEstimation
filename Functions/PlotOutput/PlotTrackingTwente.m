function [h1] = PlotTrackingTwente(datafile,nUnp,nP,Mag,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Load the data
load(datafile);
set(0, 'DefaultFigureRenderer', 'painters');

nPlotUnp = 4; % number of uneperturbed gait cycles
 


%% headers for subplots
HeadersSubplot = [];
if ~isempty(varargin)
   HeadersSubplot = varargin{1}; 
end

%% Color selection
ColGeyer = [176, 18, 52]./251;
ColCOMd = [32, 61, 179]./251;
ColGRF = [37, 186, 124]./251;
lw = 1.5;


%% Detect number of perturbation magntiudes
nPhase = length(All_GeyerCOM.Phase);
iPlots = nan(nPhase,1);
iPlots(1:nUnp) =  1;
ct = nUnp+1;
nMag = length(Mag);
for i=1:nMag
    for j = 1:nP
        iPlots(ct) = i+2 ;
        ct = ct+1;
    end
end

%% Plot the main results

h1 = figure();
set(h1,'Position',[700 181 1222 1026]);
nPhase = length(All_GeyerCOM.Phase);
t0 = 0;
dt = 0.5;
NRows = ceil(nMag/2)+1;
for i=1:nPhase
    BoolPlot = 1;
    if i==1
        t0 = 0;
    elseif iPlots(i) == iPlots(i-1)
        t0 = t0+ All_GeyerDefault.t(iSel(end)) +dt;
    else
        t0 = 0;
    end

    % select subplot
    if iPlots(i) ==1
        subplot(NRows,2,1:2);
    else
        subplot(NRows,2,iPlots(i));
    end
    if iPlots(i) ==1 && i> nPlotUnp
        BoolPlot = 0;
    end
    if BoolPlot

        % plot the tracking
        iSel = All_GeyerDefault.iDatSet(i)+1:All_GeyerDefault.iDatSet(i+1);
        plot(All_GeyerDefault.t(iSel) + t0,full(All_GeyerDefault.Tmus(iSel)),'Color',ColGeyer,'LineWidth',lw); hold on;


        % COM geyer model shooting
        %     iSel = All_GeyerCOM.iDatSet(i)+1:All_GeyerCOM.iDatSet(i+1);
        %     plot(All_GeyerCOM.t(iSel) + t0,full(All_GeyerCOM.Tmus(iSel)),'Color',ColCOMd,'LineWidth',lw);

        % COM geyer GRF model shooting
        iSel = All_GeyerCOM_GRF.iDatSet(i)+1:All_GeyerCOM_GRF.iDatSet(i+1);
        plot(All_GeyerCOM_GRF.t(iSel) + t0,full(All_GeyerCOM_GRF.Tmus(iSel)),'Color',ColGRF,'LineWidth',lw);

        % ID moment
        Set = All_GeyerCOM.Set;
        if isfield(Set,'ScaleTroque') &&Set.ScaleTroque
            Tid = All_GeyerCOM.Tid./(Set.BodyMass*Set.Height).*(70*1.75);
        else
            Tid = All_GeyerCOM.Tid;
        end
        plot(All_GeyerCOM.t(iSel) + t0, Tid(iSel),'--k');

        %     xlabel('Time [s]');
        ylabel('Moment [Nm]');
        set(gca,'FontSize',10);
        set(gca,'LineWidth',1);
        set(gca,'YLim',[-100 40]);
    end
end

subplot(NRows,2,1:2);

if ~isempty(HeadersSubplot)
    subplot(NRows,2,1:2);
    title(HeadersSubplot{1});
    for i=2:max(iPlots)-1
        subplot(NRows,2,i+1);
        title(HeadersSubplot{i});
    end
end

% add xlabel time to plot at bottom
iMax = max(iPlots);
for i=iMax-1:iMax
     subplot(NRows,2,i);
     xlabel('Time [s]');
end



% legend1 = legend({'DefaultGeyer','COMd','COMd-GRF','Tid'},'Orientation','horizontal');
legend1 = legend({'DefaultGeyer','COMd-GRF','Tid'},'Orientation','horizontal');
set(legend1,...
    'Position',[0.340837295523347 0.950531197903314 0.342880518488135 0.018518518082812],...
    'Orientation','horizontal',...
    'FontSize',10);
delete_box
end

