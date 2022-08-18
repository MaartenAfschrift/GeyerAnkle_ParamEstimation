%% Prepare data Twente for use in parameter estimation
%-----------------------------------------------------
% needed input info:

%- ankle angle
%- ankle moment
%- COM kinematics
%- COM velocity
%- Boolean to indicate gait events


% Data needed to run the scripts
% Data.RefWalk.time
% Data.RefWalk.StanceSpline_Left
% Data.RefWalk.StanceSpline_Right
% Data.RefWalk.DoubleSup_LeadLeft
% Data.RefWalk.DoubleSup_LeadRight
% Data.RefWalk.PercLeftSpline
% Data.RefWalk.PercRightSpline
% Data.RefWalk.COMdSpline
% Data.RefWalk.qSpline =  spline(tSpline,qLS');   % ankle - knee -hip
% Data.RefWalk.qRSpline =  spline(tSpline,qRS');
% Data.RefWalk.qdSpline =  spline(tSpline,qdLS');
% Data.RefWalk.qdRSpline =  spline(tSpline,qdRS');
% Data.RefWalk.Tspline =  spline(tSpline,TLS');   % ankle-knee-hip
% Data.RefWalk.TRspline =  spline(tSpline,TRS');  % ankle-knee-hip

% The main problem here is that we don't have strides in sequence. Might be
% easier to adapt the optimization code to the format we have.

% Event sequences:
%                     case 1
%                         evts = [1 2]; % HSL-TOR
%                     case 2
%                         evts = [2 3]; % TOR HSR
%                     case 3
%                         evts = [3 4]; % HSR TOL
%                     case 4
%                         evts = [4 5]; % TOL HSL
%                     case 5
%                         evts = [5 6]; % HSL TOR
%                     case 6
%                         evts = [6 7]; % TOR HSR

% info perturbation:
% perturbations were implemented as 150ms block pulses of force magnitude
% equal to 4,8,12,16% of subject body weight. perturbation were triggered
% at toe off R (based on ground reaction force).

%% load the raw datafiles
clear all; close all; clc;


Set.Freq_InterPol = 100;
Set.OutPath = 'C:\Users\u0088756\Documents\PostDoc\Data\DataTwente\DataVlutters\DataFiles_HSAdapt';
Set.ExportRefWalk = true;
Set.ExportPerturbWalk = true;
Data = load('C:\Users\u0088756\Documents\PostDoc\Data\DataTwente\DataVlutters\DataFiles_HSAdapt\SubjStatistics_all.mat');

Mass    = [ 66 89 53 62 56 80 76.3 56.6 57 66.6];
Height  = [1.79 1.92 1.62 1.82 1.80 1.94 1.89 1.73 1.78 1.70];


%%
% [nfr,~,ndir,ntr,nset,nsubj,nEvent] = size(comDatabt)

iPert   = 3;   % slow walking AP perturbations (4 = fast walking)
iDimCOM = 2;   % y values used for COM data (this seems to be AP direction)
iAnkleT = 6;   % index for selecting ankle joint moment (Left side)
iKneeT = 5;   % index for selecting ankle joint moment (Left side)
iHipT = 4;
iFoot   = 4;   % index from COM left foot
iCOMbody = 10; % index for whole body COM data
iTib    = 1;
iGas    = 2;
iEvent  = 3; % Start-HSR
tau     = 0.1;   % time delay (in s)
tauEMG  = 0.05;  % time delay reactie muscle activity
iGRFL   = 1:3;
iGRFR   = 7:9;


%% Get Data of reference walking
if Set.ExportRefWalk
    % we want all the data between heelstrike L and the subsequent heelstrike L
    % this means event sequences 1-4
    [nfr, nOut, ndim, nTR, nSet, nSubj, nEvent] = size(Data.comDataDbt);
    
    % we are only interested in Dirs 3 and 4
    
    DirNames = {'MLSLow','MLFast','Slow','Fast'};
    
    % headers of the Data structure (and splines)
    headers = {'time','COM','COMd','Ta','Tk','Th','qa','qk','qh','BoolStance','BoolSwing','DoubleSup_LeadLeft','DoubleSup_LeadRight',...
        'GRFLx','GRFLy','GRFLz','GRFRx','GRFRy','GRFRz','Tib','Gas','PertF'};
    for s = 1:nSubj
        for iPert = 3:4
            ct_Trial = 1;
            for itr=1:nTR
                % get the whole body COM data
                InputData = zeros(nfr*4,length(headers));
                t0 = 0;
                BoolErrorEvent = false;
                
                for iEvent = 1:4 % HSL-TOR-HSR-TOL-HSL
                    % com kinematics
                    COMsel = squeeze(Data.comDatabt(:,iCOMbody,iDimCOM,itr,iPert,s,iEvent));
                    COMdsel = squeeze(Data.comDataDbt(:,iCOMbody,iDimCOM,itr,iPert,s,iEvent));
                    
                    % GRF selected
                    GRFLsel = squeeze(Data.fDatabt(:,iGRFL,itr,iPert,s,iEvent));
                    GRFRsel = squeeze(Data.fDatabt(:,iGRFR,itr,iPert,s,iEvent));
                    
                    % Moments
                    Ta = squeeze(Data.jTrqDatabt(:,iAnkleT,1,itr,iPert,s,iEvent ));        % Left ankle moment
                    Tk = squeeze(Data.jTrqDatabt(:,iKneeT,1,itr,iPert,s,iEvent ));        % Left ankle moment
                    Th = squeeze(Data.jTrqDatabt(:,iHipT,1,itr,iPert,s,iEvent ));        % Left ankle moment
                    
                    % Moments are normalised
                    Ta = Ta*Mass(s)*Height(s)*9.81;
                    Tk = Tk*Mass(s)*Height(s)*9.81;
                    Th = Th*Mass(s)*Height(s)*9.81;
                    
                    % kinematics
                    qa = squeeze(Data.jAngDatabt(:,iAnkleT,1,itr,iPert,s,iEvent));
                    qk = squeeze(Data.jAngDatabt(:,iKneeT,1,itr,iPert,s,iEvent));
                    qh = squeeze(Data.jAngDatabt(:,iHipT,1,itr,iPert,s,iEvent));
                    
                    % muscle activity
                    Tib = squeeze(Data.emgDatabt(:,iTib,itr,iPert,s,iEvent));
                    Gas = squeeze(Data.emgDatabt(:,iGas,itr,iPert,s,iEvent));                                       
                    
                    % events
                    if iEvent ==1 % HSL-TOR
                        BoolStance = 1;
                        BoolSwing = 0;
                        DoubleSup_LeadLeft = 1;
                        DoubleSup_LeadRight = 0;
                    elseif iEvent == 2 % TOR HSR
                        BoolStance = 1;
                        BoolSwing = 0;
                        DoubleSup_LeadLeft = 0;
                        DoubleSup_LeadRight = 0;
                    elseif iEvent == 3 % HSR - TOL
                        BoolStance = 1;
                        BoolSwing = 0;
                        DoubleSup_LeadLeft = 0;
                        DoubleSup_LeadRight = 1;
                    else
                        BoolStance = 0;
                        BoolSwing = 1;
                        DoubleSup_LeadLeft = 0;
                        DoubleSup_LeadRight = 0;
                    end
                    
                    dt = squeeze(Data.tDatabeSeq(itr,iPert,s,iEvent));
                    if dt == 0
                        BoolErrorEvent = true;
                    end
                    iSel = (1:nfr)+(iEvent-1)*nfr;
                    InputData(iSel,1) = linspace(t0,t0+dt,nfr);
                    InputData(iSel,2) = COMsel;
                    InputData(iSel,3) = COMdsel;
                    InputData(iSel,4) = Ta;
                    InputData(iSel,5) = Tk;
                    InputData(iSel,6) = Th;
                    InputData(iSel,7) = qa;
                    InputData(iSel,8) = qk;
                    InputData(iSel,9) = qh;
                    InputData(iSel,10) = BoolStance;
                    InputData(iSel,11) = BoolSwing;
                    InputData(iSel,12) = DoubleSup_LeadLeft;
                    InputData(iSel,13) = DoubleSup_LeadRight;
                    InputData(iSel,14:16) = GRFLsel;
                    InputData(iSel,17:19) = GRFRsel;     
                    InputData(iSel,20) = Tib;    
                    InputData(iSel,21) = Gas; 
                    InputData(iSel,22) = 0;
                    t0 = t0+dt+0.000001; % small offseet for points we have twice (not very clean, but should be fine)
                end
                
                if ~BoolErrorEvent
                    % remove the NaNs (set to 0)
                    iNan = find(any(isnan(InputData)));
                    InputData(:,iNan) = 0;
                    
                    
                    % interpolate and spline the input data
                    t = InputData(:,1);
                    tSample = t(1):1/Set.Freq_InterPol:t(end);
                    DataSpline = spline(t',InputData');
                    RefWalk.(DirNames{iPert}).trial(ct_Trial).DataSpline = DataSpline;
                    RefWalk.(DirNames{iPert}).trial(ct_Trial).NaNs = iNan;
                    RefWalk.(DirNames{iPert}).trial(ct_Trial).tspan = [t(1) t(end)];
                    RefWalk.(DirNames{iPert}).trial(ct_Trial).PMag = 0;
                    RefWalk.(DirNames{iPert}).trial(ct_Trial).dtP = 0;
                    ct_Trial = ct_Trial + 1;
                end
            end
        end
        OutFolder = fullfile(Set.OutPath,['pp_' num2str(s)]);
        if ~isfolder(OutFolder)
            mkdir(OutFolder);
        end
        save(fullfile(OutFolder,'RefWalk.mat'),'RefWalk','headers');
    end
end

%% Get data of perturbed walking
PertMag = [-0.04,-0.08,-0.12,-0.16,0.04,0.08,0.12,0.16];
if Set.ExportPerturbWalk
    % we want all the data between heelstrike L and the subsequent heelstrike L
    % this means event sequences 1-4
    [nfr, nOut, ndim, nTR, nP, nSet, nSubj, nEvent] = size(Data.comDataDpt);
    
    % we are only interested in Dirs 3 and 4
    
    DirNames = {'MLSLow','MLFast','Slow','Fast'};
    
    % headers of the Data structure (and splines)
    headers = {'time','COM','COMd','Ta','Tk','Th','qa','qk','qh','BoolStance','BoolSwing','DoubleSup_LeadLeft','DoubleSup_LeadRight',...
        'GRFLx','GRFLy','GRFLz','GRFRx','GRFRy','GRFRz','Tib','Gas','PertF'};
    for s = 1:nSubj     % loop subjects
        for iPert = 3:4     % loop walkings speeds
            for iP = 1:nP       % loop perturabtion profiles
                ct_Trial = 1;
                for itr=1:nTR       % number of trials
                    % get the whole body COM data
                    InputData = zeros(nfr*4,13);
                    t0 = 0;
                    BoolErrorEvent = false;
                    
                    for iEvent = 1:4    % 4 event sequences (i.e. two consequtive heelstrikes)
                        % com kinematics
                        COMsel = squeeze(Data.comDatapt(:,iCOMbody,iDimCOM,itr,iP,iPert,s,iEvent));
                        COMdsel = squeeze(Data.comDataDpt(:,iCOMbody,iDimCOM,itr,iP,iPert,s,iEvent));
                        
                        % GRF selected
                        GRFLsel = squeeze(Data.fDatapt(:,iGRFL,itr,iP,iPert,s,iEvent));
                        GRFRsel = squeeze(Data.fDatapt(:,iGRFR,itr,iP,iPert,s,iEvent));
                        
                        % Moments
                        Ta = squeeze(Data.jTrqDatapt(:,iAnkleT,1,itr,iP,iPert,s,iEvent ));        % Left ankle moment
                        Tk = squeeze(Data.jTrqDatapt(:,iKneeT,1,itr,iP,iPert,s,iEvent ));        % Left ankle moment
                        Th = squeeze(Data.jTrqDatapt(:,iHipT,1,itr,iP,iPert,s,iEvent ));        % Left ankle moment
                        
                        % Moments are normalised
                        Ta = Ta*Mass(s)*Height(s)*9.81;
                        Tk = Tk*Mass(s)*Height(s)*9.81;
                        Th = Th*Mass(s)*Height(s)*9.81;
                        
                        % kinematics
                        qa = squeeze(Data.jAngDatapt(:,iAnkleT,1,itr,iP,iPert,s,iEvent));
                        qk = squeeze(Data.jAngDatapt(:,iKneeT,1,itr,iP,iPert,s,iEvent));
                        qh = squeeze(Data.jAngDatapt(:,iHipT,1,itr,iP,iPert,s,iEvent));
                        
                        % muscle activity
                        Tib = squeeze(Data.emgDatapt(:,iTib,itr,iP,iPert,s,iEvent));
                        Gas = squeeze(Data.emgDatapt(:,iGas,itr,iP,iPert,s,iEvent));
                        
                        % events
                        if iEvent ==1 % HSL-TOR
                            BoolStance = 1;
                            BoolSwing = 0;
                            DoubleSup_LeadLeft = 1;
                            DoubleSup_LeadRight = 0;
                        elseif iEvent == 2 % TOR HSR
                            BoolStance = 1;
                            BoolSwing = 0;
                            DoubleSup_LeadLeft = 0;
                            DoubleSup_LeadRight = 0;
                        elseif iEvent == 3 % HSR - TOL
                            BoolStance = 1;
                            BoolSwing = 0;
                            DoubleSup_LeadLeft = 0;
                            DoubleSup_LeadRight = 1;
                        else
                            BoolStance = 0;
                            BoolSwing = 1;
                            DoubleSup_LeadLeft = 0;
                            DoubleSup_LeadRight = 0;
                        end
                        
                        dt = squeeze(Data.tDatapeSeq(itr,iP,iPert,s,iEvent));
                        if dt == 0
                            BoolErrorEvent = true;
                        end
                        
                        % perturbation force
                        PForce = zeros(length(Ta),1);
                        if iEvent == 2
                            % select indices with perturbation (first 150ms
                            % after event)
                            t = linspace(t0,t0+dt,nfr);
                            iPertSel = find(t-t(1)<=0.150);
                            PForce(iPertSel) = PertMag(iP).*Mass(s).*9.81;
                            if ~isempty(iPertSel)                            
                                t0P = t(iPertSel(1));
                                tendP = t(iPertSel(end));
                                dtP = [t0P tendP];
                            else
                                dtP = NaN;
                            end
                        end
                        
                        iSel = (1:nfr)+(iEvent-1)*nfr;
                        InputData(iSel,1) = linspace(t0,t0+dt,nfr);
                        InputData(iSel,2) = COMsel;
                        InputData(iSel,3) = COMdsel;
                        InputData(iSel,4) = Ta;
                        InputData(iSel,5) = Tk;
                        InputData(iSel,6) = Th;
                        InputData(iSel,7) = qa;
                        InputData(iSel,8) = qk;
                        InputData(iSel,9) = qh;
                        InputData(iSel,10) = BoolStance;
                        InputData(iSel,11) = BoolSwing;
                        InputData(iSel,12) = DoubleSup_LeadLeft;
                        InputData(iSel,13) = DoubleSup_LeadRight;
                        InputData(iSel,14:16) =  GRFLsel;
                        InputData(iSel,17:19) =  GRFRsel;
                        InputData(iSel,20) =  Tib;
                        InputData(iSel,21) =  Gas;
                        InputData(iSel,22) =  PForce;
                        t0 = t0+dt+0.000001; % small offseet for points we have twice (not very clean, but should be fine)
                    end
                    
                    if ~BoolErrorEvent && ~any(isnan(InputData(:,1)))
                        % remove the NaNs (set to 0)
                        iNan = find(any(isnan(InputData)));
                        InputData(:,iNan) = 0;
                        
                        
                        % interpolate and spline the input data
                        t = InputData(:,1);
                        tSample = t(1):1/Set.Freq_InterPol:t(end);
                        DataSpline = spline(t',InputData');
                        PertWalk.(DirNames{iPert}).PType(iP).trial(ct_Trial).DataSpline = DataSpline;
                        PertWalk.(DirNames{iPert}).PType(iP).trial(ct_Trial).NaNs = iNan;
                        PertWalk.(DirNames{iPert}).PType(iP).trial(ct_Trial).tspan = [t(1) t(end)];
                        PertWalk.(DirNames{iPert}).PType(iP).trial(ct_Trial).PMag = PertMag(iP).*Mass(s).*9.81;
                        PertWalk.(DirNames{iPert}).PType(iP).trial(ct_Trial).dtP = dtP;
                        ct_Trial = ct_Trial + 1;
                    else
                        disp('skipped trial');
                    end
                end
            end
        end
        OutFolder = fullfile(Set.OutPath,['pp_' num2str(s)]);
        if ~isfolder(OutFolder)
            mkdir(OutFolder);
        end
        save(fullfile(OutFolder,'PertWalk.mat'),'PertWalk','headers');
    end
end