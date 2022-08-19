function [h1] = PlotTrackingTwente_MuscleAv_Exp(datafile,RefWalk,PertWalk,nUnp,nP,Mag,varargin)
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

%% get the measured muscle activity

% reconstruct the Datastructure used for data analysis
Set = All_GeyerCOM.Set;
% unperturbed walking: Get at least n cycles without missing data
i = 16;
CombDatSel = [];
ct = 1;

iTrials = TrialInfo(FileInfo == 0);
CombDatSel = RefWalk.(Set.WalkSpeed).trial(iTrials );
for iP = Set.iPertSel
    iTrials = TrialInfo(FileInfo == iP);
    DatP = PertWalk.(Set.WalkSpeed).PType(iP).trial(iTrials);
    CombDatSel = [CombDatSel DatP];
end

%% Compute the average muscle activity for each perturbation type

h1 = figure();
set(h1,'Position',[462      423     1135    463]);   
nPhase = length(All_GeyerCOM.Phase);
t0 = 0;
dt = 0.5;
NRows = ceil(nMag/2)+1;

AllMusAct = nan(100,nPhase,3,2);
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
    act_Geyer = full(All_GeyerDefault.x(1:2,iSel)) ;   
    AllMusAct(:,i,1,:) = interp1(1:length(act_Geyer),act_Geyer',linspace(0,length(act_Geyer),100));
        
    % COM geyer GRF model shooting
    iSel = All_GeyerCOM_GRF.iDatSet(i)+1:All_GeyerCOM_GRF.iDatSet(i+1);
    act_GRF = full(All_GeyerCOM_GRF.x(1:2,iSel));
    AllMusAct(:,i,2,:) = interp1(1:length(act_GRF),act_GRF',linspace(0,length(act_GRF),100)); 
    
    % get the experimental data
    tVect = linspace(CombDatSel(i).tspan(1),CombDatSel(i).tspan(2),100);
    dsel = ppval(tVect,CombDatSel(i).DataSpline)';
    AllMusAct(:,i,3,:) = dsel(:,[21 20]);
    
end

% compute averages
iP_Unique = unique(iPlots);
MuscleAct = nan(100,length(iP_Unique),3,2); 
for i=1:length(iP_Unique)
    % find indices for this perturbation
    iS = find(iPlots==iP_Unique(i));
    % compute average
    MuscleAct(:,i,:,:) =  nanmean(AllMusAct(:,iS,:,:),2);      
end
ColSel = [];

CRef = [0.4 0.4 0.4];
CVectBack    =  [182 221 255; 113 188 253;  38 129 208;    27 90 145]./255;
CVectFor   =  [255 214 165; 247 188 127;  189 50 60;     129 32 38]./255;
Cols = [CRef;  CVectBack; CVectFor ];

for i=1:3    
    for m=1:2
        subplot(2,3,(m-1)*3+i);
        plot(MuscleAct(:,1,i,m),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
        for j=2:length(iP_Unique)
            plot(MuscleAct(:,j,i,m),'Color',Cols(j,:),'LineWidth',lw); hold on;
        end
        plot(MuscleAct(:,1,i,m),'Color',Cols(1,:),'LineWidth',lw*2); hold on;
        if i<3
            set(gca,'YLim',[0 0.5]);
            set(gca,'YTick',[0 0.25 0.5]);
        else
            
        end
        set(gca,'XTick',[0 50 100]);
        set(gca,'FontSize',10);
        set(gca,'LineWidth',2);
    end
end
legend1 = legend(HeadersSubplot);
set(legend1,...
    'Position',[0.854797549118755 0.250491949488068 0.108370042162319 0.269978394196821]);

subplot(2,3,1)
ylabel('Soleus activity []');
title('Default Geyer');
subplot(2,3,2)
title('COMd-GRF');
subplot(2,3,3)
title('EMG ');
ylabel('Gastrocn. EMG []');

subplot(2,3,4)
ylabel('Tibialis activity []');
xlabel('% gait cycle');
subplot(2,3,5);
xlabel('% gait cycle');
subplot(2,3,6);
ylabel('Tibialis EMG []');
xlabel('% gait cycle');


delete_box


end

