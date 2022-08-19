function [h1] = PlotTrackingTwente_Soleus(datafile,nUnp,nP,Mag,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Load the data
load(datafile);

%% headers for subplots
HeadersSubplot = [];
if ~isempty(varargin)
   HeadersSubplot = varargin{1}; 
end

%% Color selection
ColGeyer = [37, 186, 124]./251;
ColCOMd = [32, 61, 179]./251;
ColGRF = [176, 18, 52]./251;
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
    
    % plot the Default tracking
    iSel = All_GeyerDefault.iDatSet(i)+1:All_GeyerDefault.iDatSet(i+1);
    plot(All_GeyerDefault.t(iSel) + t0,full(All_GeyerDefault.x(1,iSel)),'Color',ColGeyer,'LineWidth',lw); hold on;
        
    % COM geyer GRF model shooting
    iSel = All_GeyerCOM_GRF.iDatSet(i)+1:All_GeyerCOM_GRF.iDatSet(i+1);
    plot(All_GeyerCOM_GRF.t(iSel) + t0,full(All_GeyerCOM_GRF.x(1,iSel)),'Color',ColGRF,'LineWidth',lw);
    
    % ID moment
    %plot(All_GeyerCOM.t(iSel) + t0, All_GeyerCOM.Tid(iSel),'--k');   
    
%     xlabel('Time [s]');
    ylabel('Moment [Nm]');
    set(gca,'FontSize',10);
    set(gca,'LineWidth',1);
    set(gca,'YLim',[0 1]);
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



legend1 = legend({'DefaultGeyer','COMd-GRF'},'Orientation','horizontal');
set(legend1,...
    'Position',[0.340837295523347 0.950531197903314 0.342880518488135 0.018518518082812],...
    'Orientation','horizontal',...
    'FontSize',10);
delete_box
end

