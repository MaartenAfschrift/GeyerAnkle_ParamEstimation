function [h1] = PlotTrackingTwente_TorqueAv_Exo(datafile,nUnp,nP,Mag,varargin)
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

%% Compute the average muscle activity for each perturbation type

h1 = figure();
set(h1,'Position',[176      380     1466    294]);
nPhase = length(All_GeyerCOM.Phase);
t0 = 0;
dt = 0.5;
NRows = ceil(nMag/2)+1;

AllDatVect = nan(100,nPhase,3);
for i=1:nPhase
    
    if i==1
        t0 = 0;
    elseif iPlots(i) == iPlots(i-1)
        t0 = t0+ All_GeyerDefault.t(iSel(end)) +dt;
    else
        t0 = 0;
    end
    
    % plot the Default tracking
    iSel = All_GeyerDefault.iDatSet(i)+1:All_GeyerDefault.iDatSet(i+1);
    act_Geyer = full(All_GeyerDefault.Tmus(iSel)) ;
    AllDatVect(:,i,1) = interp1(1:length(act_Geyer),act_Geyer',linspace(0,length(act_Geyer),100));
    
    % COM geyer GRF model shooting
    iSel = All_GeyerCOM_GRF.iDatSet(i)+1:All_GeyerCOM_GRF.iDatSet(i+1);
    act_GRF = full(All_GeyerCOM_GRF.Tmus(iSel));
    AllDatVect(:,i,2) = interp1(1:length(act_GRF),act_GRF',linspace(0,length(act_GRF),100));
    
    Set = All_GeyerCOM_GRF.Set;
    if isfield(Set,'ScaleTroque') &&Set.ScaleTroque    
        Tid = All_GeyerCOM_GRF.Tid./(Set.BodyMass*Set.Height).*(70*1.75);
    else
        Tid = All_GeyerCOM_GRF.Tid;
    end
    
    act_Tid = full(Tid(iSel));
    AllDatVect(:,i,3) = interp1(1:length(act_Tid),act_Tid',linspace(0,length(act_Tid),100));
end

% compute averages
iP_Unique = unique(iPlots);
AllDat = nan(100,length(iP_Unique),3);
ExoDat = nan(100,length(iP_Unique));
for i=1:length(iP_Unique)
    % find indices for this perturbation
    iS = find(iPlots==iP_Unique(i));
    % compute average
    AllDat(:,i,:) =  nanmean(AllDatVect(:,iS,:),2);
    ContDiff = AllDatVect(:,iS,2) - AllDatVect(:,iS,1);
    ExoDat(:,i) = nanmean(ContDiff,2);
end



%% Plot the figure
CRef = [0.4 0.4 0.4];
CVectBack    =  [182 221 255; 113 188 253;  38 129 208;    27 90 145]./255;
CVectFor   =  [255 214 165; 247 188 127;  189 50 60;     129 32 38]./255;
Cols = [CRef;  CVectBack; CVectFor ];

% experimental alues
subplot(1,3,1);
plot(AllDat(:,1,3),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
for j=2:length(iP_Unique)
    plot(AllDat(:,j,3),'Color',Cols(j,:),'LineWidth',lw); hold on;
end
plot(AllDat(:,1,3),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
set(gca,'YLim',[-100 50]);
set(gca,'YTick',[-100 -50 0 50]);

% delta torque experiment
subplot(1,3,2);
plot(AllDat(:,1,3) - AllDat(:,1,3),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
for j=2:length(iP_Unique)
    plot(AllDat(:,j,3) - AllDat(:,1,3),'Color',Cols(j,:),'LineWidth',lw); hold on;
end
set(gca,'YLim',[-100 100]);
set(gca,'YTick',[-100 -50 0 50 100]);

% exoskeleton controller
subplot(1,3,3);
plot(ExoDat(:,1),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
for j=2:length(iP_Unique)
    plot(ExoDat(:,j),'Color',Cols(j,:),'LineWidth',lw); hold on;
end
set(gca,'YLim',[-100 100]);
set(gca,'YTick',[-100 -50 0 50 100]);

for i=1:3
    subplot(1,3,i)
    set(gca,'XTick',[0 50 100]);
    set(gca,'FontSize',10);
    set(gca,'LineWidth',2);
end

legend1 = legend(HeadersSubplot);
set(legend1,...
    'Position',[0.912512469206906 0.448499685434533 0.0839017720697356 0.425170056167102]);

subplot(1,3,1)
ylabel('Ankle Moment [Nm]');
title('Experiment');

subplot(1,3,2)
ylabel('\Delta Ankle Moment [Nm]');
title('Experiment');

subplot(1,3,3)
ylabel('Desired Exo Torque [Nm]');
title('Proposed Controller');

for i=1:3
    xlabel('% gait cycle');
end

delete_box


end

