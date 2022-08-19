%% Create casadi function for integration in advance

% this function creates a bunch of casadifunctions for the forward simulation. 
% the idea is the create the symbolic equations here for 1000 fixed integration steps
% these casadi expression are saved and loaded in the problem formulation.

% the only reason why we do this here is the speed-up the problem formulation in 
% GainEstimation_AnkleOnly_nCycl_Shooting.m

% clear all variables
clear all; close all; clc;

%% settings
% path to save the casadifiles
MainPath    = 'C:\Users\Maarten\Documents\Software\Sim\GeyerAnkle_ParamEstimation';
OutPath     = fullfile(MainPath,'Functions','CasadiFuncs');

Sets =  {'GeyerDefault','GeyerDefault_Limited1','Geyer_COMd',...
'Geyer_COMd_Limited2','Geyer_COMd_Limited1','Geyer_COMd_GRF',...
'Geyer_Default_GRF','Geyer_COMd_GRF_lim','Geyer_COMdd_GRF'};
% Sets =  {'GeyerDefault'};
% Sets = {'Geyer_Default_GRF'};
% Sets = {'Geyer_COMd_GRF_lim'};
% Sets = {'Geyer_COMAll_GRF','Geyer_COMdd_GRF'};
for iSet=1:length(Sets)

    Set.GainsEstimation  = Sets{iSet};
    Set.ATendon = [20; 35];
    
    % general muscle information
    nMuscles    = 2;
    iSol        = 1;
    iTib        = 2;
    iMuscle     = [8 9];
    ATendon     = Set.ATendon;
    shift       = getShift(ATendon);
    
    % path information    
    ModelPath       = fullfile(MainPath,'osimModel','Models','gait_2D_Strength.osim');
    
    % get the muscle parameters
    MuscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
        'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
    [params]=ReadMuscleParameters(ModelPath,MuscleNames);
    params = params(:,iMuscle);
    
    %% Create function to compute state derivative at gait phases
    
    % casadi libraries
    import casadi.*
    
    % initial state
    x = SX.sym('x',6,1);
    
    % feedback gains (depending on the parameter estimation case used)
    if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        OptReflex    = SX.sym('OptReflex',6);
        Sol.e0      = OptReflex(1);
        Sol.G       = OptReflex(2);
        Tib.e0      = OptReflex(3);
        Tib.G       = OptReflex(4);
        Tib.loff    = OptReflex(5);
        Tib.GSol    = OptReflex(6);
        Sol.COMd    = 0;
        Tib.COMd    = 0;
    elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
        OptReflex   = SX.sym('OptReflex',10);
        Sol.e0      = OptReflex(1);
        Sol.G       = OptReflex(2);
        Tib.e0      = OptReflex(3);
        Tib.G       = OptReflex(4);
        Tib.loff    = OptReflex(5);
        Tib.GSol    = OptReflex(6);
        Sol.COM     = OptReflex(7);
        Sol.COMd    = OptReflex(8);
        Tib.COM     = OptReflex(9);
        Tib.COMd    = OptReflex(10);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
        OptReflex    = SX.sym('OptReflex',8);
        Sol.e0      = OptReflex(1);
        Sol.G       = OptReflex(2);
        Tib.e0      = OptReflex(3);
        Tib.G       = OptReflex(4);
        Tib.loff    = OptReflex(5);
        Tib.GSol    = OptReflex(6);
        Sol.COMd    = OptReflex(7);
        Tib.COMd    = OptReflex(8);
    elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')
        OptReflex   = SX.sym('OptReflex',3);
        ConstReflex = SX.sym('ConstReflex',3);
        Sol.e0      = ConstReflex(1);
        Sol.G       = OptReflex(1);
        Tib.e0      = ConstReflex(2);
        Tib.G       = OptReflex(2);
        Tib.loff    = ConstReflex(3);
        Tib.GSol    = OptReflex(3);
        Sol.COMd    = 0;
        Tib.COMd    = 0;
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
        OptReflex    = SX.sym('OptReflex',5);
        ConstReflex = SX.sym('ConstReflex',3);
        Sol.e0      = ConstReflex(1);
        Sol.G       = OptReflex(1);
        Tib.e0      = ConstReflex(2);
        Tib.G       = OptReflex(2);
        Tib.loff    = ConstReflex(3);
        Tib.GSol    = OptReflex(3);
        Sol.COMd    = OptReflex(4);
        Tib.COMd    = OptReflex(5);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        OptReflex    = SX.sym('OptReflex',6);
        ConstReflex = SX.sym('ConstReflex',2);
        Sol.e0      = ConstReflex(1);
        Sol.G       = OptReflex(1);
        Tib.e0      = OptReflex(2);
        Tib.G       = OptReflex(3);
        Tib.loff    = ConstReflex(2);
        Tib.GSol    = OptReflex(4);
        Sol.COMd    = OptReflex(5);
        Tib.COMd    = OptReflex(6);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        OptReflex    = SX.sym('OptReflex',7);
        ConstReflex = SX.sym('ConstReflex',1);
        Sol.e0      = ConstReflex(1);
        Sol.G       = OptReflex(1);
        Tib.e0      = OptReflex(2);
        Tib.G       = OptReflex(3);
        Tib.loff    = OptReflex(4);
        Tib.GSol    = OptReflex(5);
        Sol.COMd    = OptReflex(6);
        Tib.COMd    = OptReflex(7);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        OptReflex    = SX.sym('OptReflex',10);
        Sol.e0      = OptReflex(1);
        Sol.G       = OptReflex(2);
        Tib.e0      = OptReflex(3);
        Tib.G       = OptReflex(4);
        Tib.loff    = OptReflex(5);
        Tib.GSol    = OptReflex(6);
        Sol.COMd    = OptReflex(7);
        Tib.COMd    = OptReflex(8);
        Sol.COMdd   = OptReflex(9);
        Tib.COMdd   = OptReflex(10);
    elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        OptReflex    = SX.sym('OptReflex',12);
        Sol.e0      = OptReflex(1);
        Sol.G       = OptReflex(2);
        Tib.e0      = OptReflex(3);
        Tib.G       = OptReflex(4);
        Tib.loff    = OptReflex(5);
        Tib.GSol    = OptReflex(6);
        Sol.COM     = OptReflex(7);
        Tib.COM     = OptReflex(8);
        Sol.COMd    = OptReflex(9);
        Tib.COMd    = OptReflex(10);
        Sol.COMdd   = OptReflex(11);
        Tib.COMdd   = OptReflex(12);
    end
    
    % Muscle tendon length and COM velocity as input arguments in this function
    lMT = SX.sym('lMT',nMuscles,1);
    COMd_delay = SX.sym('COMd_delay',1,1);
    
    % optimal input arguments (COM_delay and COMdd_delay)
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        COMdd_delay = SX.sym('COMdd_delay',1,1);
    end
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') 
        COM_delay = SX.sym('COM_delay',1,1);
    end

    % GRFy if this is used to normalise soleus reflex
    if strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF') ||...
            strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim') || strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        GRFy = SX.sym('GRFy',1,1);
    end
    
    % unpack the state vector
    a = x(1:2);
    lMtilde = x(3:4);
    yPade_Fb = x(5);
    yPade_lM = x(6);
    
    % explicit contraction dynamics
    [dlMtildedt,FT] = Explicit_ContractionDyn_limitFMV(a,lMtilde,lMT,params,ATendon,shift);
    
    % pade aproximation
    Fmus_norm = FT(1)./params(1,1);
    [Fmus_t,yPade_Fb_dot] = PadeApprox_order1_20ms(Fmus_norm,yPade_Fb);
    [lM_t,yPade_lM_dot] = PadeApprox_order1_20ms(lMtilde(iTib),yPade_lM);
        
    % --------------------------------
    % --------- Stance phase----------
    %---------------------------------    
    % stance phase reflexes
    % Soleus Reflex
    if strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim') || strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        eSol = Sol.G.*GRFy.*Fmus_t;
    else
        eSol = Sol.G.*Fmus_t;
    end
    
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        COMfb = COM_delay.*Sol.COM + COMd_delay.*Sol.COMd + COMdd_delay.*Sol.COMdd;
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')     
        COMfb = COMd_delay.*Sol.COMd + COMdd_delay.*Sol.COMdd;
    else
        COMfb = COMd_delay.*Sol.COMd;
    end
    eSol = eSol + COMfb;
    
    % Tibialis anterior reflex
    eTib   = Tib.G.*(lM_t-Tib.loff) - Tib.GSol.*Fmus_t;
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        COMfb = COM_delay.*Tib.COM + COMd_delay.*Tib.COMd + COMdd_delay.*Tib.COMdd;
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')     
        COMfb = COMd_delay.*Tib.COMd + COMdd_delay.*Tib.COMdd;
    else
        COMfb = COMd_delay.* Tib.COMd;
    end
    eTib = eTib + COMfb;
    
    % combine excitations in one vector
    eFB = [eSol; eTib];
    % allow saturation of values between 0 and 1
    [eFB] = SaturateExcitation(eFB,15);
    % add baseline activity after saturation
    eBaseline = [Sol.e0 Tib.e0]';
    eFB = eFB + eBaseline;
    % muscle activation dynamics
    dadt = ActivationDynamics(eFB,a,0.02,0.04,0.1);
    % get state derivative
    xdot = [ dadt; dlMtildedt; yPade_Fb_dot; yPade_lM_dot];
    % create the casadi function
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,ConstReflex,x,lMT,COMd_delay},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,x,lMT,COMd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,ConstReflex,x,lMT,COMd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,x,lMT,COM_delay,COMd_delay,COMdd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,x,lMT,COMd_delay,COMdd_delay,GRFy},{xdot,FT});
    else
        f_Xdot_Stance = Function('f_Xdot_Stance',{OptReflex,x,lMT,COMd_delay},{xdot,FT});
    end
    
    
    % --------------------------------
    % --------- Double Support -------
    %---------------------------------
    % stance phase reflexes
    % Soleus Reflex
    if strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim') || strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        eSol = Sol.G.*GRFy.*Fmus_t;
    else
        eSol = Sol.G.*Fmus_t;
    end
    
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        COMfb = COM_delay.*Sol.COM + COMd_delay.*Sol.COMd + COMdd_delay.*Sol.COMdd;
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')     
        COMfb = COMd_delay.*Sol.COMd + COMdd_delay.*Sol.COMdd;
    else
        COMfb = COMd_delay.*Sol.COMd;
    end
    eSol = eSol + COMfb;
    
    % Tibialis anterior reflex
    eTib   = Tib.G.*(lM_t-Tib.loff) - Tib.GSol.*Fmus_t;
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        COMfb = COM_delay.*Tib.COM + COMd_delay.*Tib.COMd + COMdd_delay.*Tib.COMdd;
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')     
        COMfb = COMd_delay.*Tib.COMd + COMdd_delay.*Tib.COMdd;
    else
        COMfb = COMd_delay.* Tib.COMd;
    end
    eTib = eTib + COMfb;
    % combine excitations in one vector
    eFB = [eSol; eTib];
    % allow saturation of values between 0 and 1
    [eFB] = SaturateExcitation(eFB,15);
    % add baseline activity after saturation
    eBaseline = [Sol.e0 Tib.e0]';
    eFB = eFB + eBaseline;
    % muscle activation dynamics
    dadt = ActivationDynamics(eFB,a,0.02,0.04,0.1);
    % get state derivative
    xdot = [ dadt; dlMtildedt; yPade_Fb_dot; yPade_lM_dot];
    % create the casadi function
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,ConstReflex,x,lMT,COMd_delay},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,x,lMT,COMd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,ConstReflex,x,lMT,COMd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,x,lMT,COM_delay,COMd_delay,COMdd_delay,GRFy},{xdot,FT});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,x,lMT,COMd_delay,COMdd_delay,GRFy},{xdot,FT});
    else
        f_Xdot_DS = Function('f_Xdot_DS',{OptReflex,x,lMT,COMd_delay},{xdot,FT});
    end
    
    % --------------------------------
    % --------- Swing Phase--- -------
    %---------------------------------
    
    eSol = 0;
    % Tibialis anterior reflex
    eTib   = Tib.G.*(lM_t-Tib.loff);
    % combine excitations in one vector
    eFB = [eSol; eTib];
    % allow saturation of values between 0 and 1
    [eFB] = SaturateExcitation(eFB,15);
    % add baseline activity after saturation
    eBaseline = [Sol.e0 Tib.e0]';
    eFB = eFB + eBaseline;
    % muscle activation dynamics
    dadt = ActivationDynamics(eFB,a,0.02,0.04,0.1);
    % get state derivative
    xdot = [ dadt; dlMtildedt; yPade_Fb_dot; yPade_lM_dot];
    % create the casadi function
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Xdot_Swing = Function('f_Xdot_Swing',{OptReflex,ConstReflex,x,lMT},{xdot,FT});
    else
        f_Xdot_Swing = Function('f_Xdot_Swing',{OptReflex,x,lMT},{xdot,FT});
    end
       
    
    %%      Create Integration scheme -- Stance Phase
    %--------------------------------------------------
    %--------------------------------------------------
    
    Ntot                = 1000;
    h                   = 1./Ntot;
    States              = SX(6,Ntot);
    J                   = SX(Ntot,1);
    TMusV               = SX(Ntot,1);
    lMT                 = SX.sym('lMT',2,Ntot);
    COMd_delay          = SX.sym('COMd_delay',1,Ntot);
    dMs                 = SX.sym('dMs',2,Ntot);
    x0                  = SX.sym('x',6,1);
    Tid                 = SX.sym('Tid',1,Ntot);
    if strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim') || strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        GRFy            = SX.sym('lMT',1,Ntot);
    end
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        COM_delay      = SX.sym('COM_delay',1,Ntot);
    end
    if strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF') || strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF') 
        COMdd_delay      = SX.sym('COMdd_delay',1,Ntot);
    end
    
    States(:,1)         = x0;
    ctx = 1;
    for i = 1:Ntot-1
        % unpack the state vector
        a = States(1:2,ctx);
        lMtilde = States(3:4,ctx);
        yPade_Fb = States(5,ctx);
        yPade_lM = States(6,ctx);
        
        x = [a; lMtilde; yPade_Fb; yPade_lM];
        
        % compute the state derivative
        if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || ...
                strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
            [xdot,FT] = f_Xdot_Stance(OptReflex,ConstReflex,x,lMT(:,i),COMd_delay(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
            [xdot,FT] = f_Xdot_Stance(OptReflex,x,lMT(:,i),COMd_delay(i),GRFy(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
            [xdot,FT] = f_Xdot_Stance(OptReflex,ConstReflex,x,lMT(:,i),COMd_delay(i),GRFy(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
            [xdot,FT] = f_Xdot_Stance(OptReflex,x,lMT(:,i),COM_delay(i),COMd_delay(i),COMdd_delay(i),GRFy(i)); 
        elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
            [xdot,FT] = f_Xdot_Stance(OptReflex,x,lMT(:,i),COMd_delay(i),COMdd_delay(i),GRFy(i)); 
        else
            [xdot,FT] = f_Xdot_Stance(OptReflex,x,lMT(:,i),COMd_delay(i));
        end
        
        dadt = xdot(1:2);
        dlMtildedt = xdot(3:4);
        yPade_Fb_dot = xdot(5);
        yPade_lM_dot = xdot(6);
        
        % forward euler integrator
        x_step = [a + dadt*h; lMtilde + dlMtildedt*h; yPade_Fb + yPade_Fb_dot*h; yPade_lM + yPade_lM_dot*h];
        States(:,ctx+1) = x_step;
        % evaluate objective function
        TMus = sum(dMs(:,i).*FT);
        Terror = Tid(i) - TMus;
        J(ctx) = (Terror/10)^2;
        
        % aux variables
        TMusV(ctx) = TMus;
        
        % append counter states
        ctx = ctx+1;
    end
    
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,x0,lMT,COM_delay,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,x0,lMT,COM_delay,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,x0,lMT,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,x0,lMT,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{J});
    else
        f_Int1000_Stance_Full = Function('f_Int1000_Stance_Full',{OptReflex,x0,lMT,COMd_delay,dMs,Tid},{States,J,TMusV});
        f_Int1000_Stance_J = Function('f_Int1000_Stance_J',{OptReflex,x0,lMT,COMd_delay,dMs,Tid},{J});
    end
    if ~isfolder(OutPath)
        mkdir(OutPath);
    end
    f_Int1000_Stance_Full.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_Stance_Full']));
    f_Int1000_Stance_J.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_Stance_J']));
    
    clear States J TmusV

    %----------------------------------------------
    %%  Create Integration scheme -- Double support
    %----------------------------------------------
    %----------------------------------------------
    
    
    States              = SX(6,Ntot);
    J                   = SX(Ntot,1);
    TMusV               = SX(Ntot,1);
    States(:,1)         = x0;
    ctx = 1;
    for i = 1:Ntot-1
        % unpack the state vector
        a = States(1:2,ctx);
        lMtilde = States(3:4,ctx);
        yPade_Fb = States(5,ctx);
        yPade_lM = States(6,ctx);
        
        x = [a; lMtilde; yPade_Fb; yPade_lM];
        
        % compute the state derivative
        if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || ...
                strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
            [xdot,FT] = f_Xdot_DS(OptReflex,ConstReflex,x,lMT(:,i),COMd_delay(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
            [xdot,FT] = f_Xdot_DS(OptReflex,x,lMT(:,i),COMd_delay(i),GRFy(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
            [xdot,FT] = f_Xdot_DS(OptReflex,ConstReflex,x,lMT(:,i),COMd_delay(i),GRFy(i));
        elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
            [xdot,FT] = f_Xdot_DS(OptReflex,x,lMT(:,i),COM_delay(i),COMd_delay(i),COMdd_delay(i),GRFy(i)); 
        elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
            [xdot,FT] = f_Xdot_DS(OptReflex,x,lMT(:,i),COMd_delay(i),COMdd_delay(i),GRFy(i)); 
        else
            [xdot,FT] = f_Xdot_DS(OptReflex,x,lMT(:,i),COMd_delay(i));
        end
            
        
        dadt = xdot(1:2);
        dlMtildedt = xdot(3:4);
        yPade_Fb_dot = xdot(5);
        yPade_lM_dot = xdot(6);
        
        % forward euler integrator
        x_step = [a + dadt*h; lMtilde + dlMtildedt*h; yPade_Fb + yPade_Fb_dot*h; yPade_lM + yPade_lM_dot*h];
        States(:,ctx+1) = x_step;
        % evaluate objective function
        TMus = sum(dMs(:,i).*FT);
        Terror = Tid(i) - TMus;
        J(ctx) = (Terror/10)^2;
        
        % aux variables
        TMusV(ctx) = TMus;
        
        % append counter states
        ctx = ctx+1;
    end
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid},{J});
   elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,ConstReflex,x0,lMT,COMd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,x0,lMT,COM_delay,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,x0,lMT,COM_delay,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{J});
    elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,x0,lMT,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,x0,lMT,COMd_delay,COMdd_delay,dMs,Tid,GRFy},{J});
    else
        f_Int1000_DS_Full = Function('f_Int1000_DS_Full',{OptReflex,x0,lMT,COMd_delay,dMs,Tid},{States,J,TMusV});
        f_Int1000_DS_J = Function('f_Int1000_DS_J',{OptReflex,x0,lMT,COMd_delay,dMs,Tid},{J});
    end
    f_Int1000_DS_Full.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_DS_Full']));
    f_Int1000_DS_J.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_DS_J']));
    clear States J TmusV

    %----------------------------------------------
    %%  Create Integration scheme -- Swing Phase
    %----------------------------------------------
    %----------------------------------------------
    
    
    States              = SX(6,Ntot);
    J                   = SX(Ntot,1);
    TMusV               = SX(Ntot,1);
    States(:,1)         = x0;
    ctx = 1;
    for i = 1:Ntot-1
        % unpack the state vector
        a = States(1:2,ctx);
        lMtilde = States(3:4,ctx);
        yPade_Fb = States(5,ctx);
        yPade_lM = States(6,ctx);
        
        x = [a; lMtilde; yPade_Fb; yPade_lM];
        
        % compute the state derivative
        if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') ||...
                strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
            [xdot,FT] =f_Xdot_Swing(OptReflex,ConstReflex,x,lMT(:,i));
        else
            [xdot,FT] =f_Xdot_Swing(OptReflex,x,lMT(:,i));
        end
        
        dadt = xdot(1:2);
        dlMtildedt = xdot(3:4);
        yPade_Fb_dot = xdot(5);
        yPade_lM_dot = xdot(6);
        
        % forward euler integrator
        x_step = [a + dadt*h; lMtilde + dlMtildedt*h; yPade_Fb + yPade_Fb_dot*h; yPade_lM + yPade_lM_dot*h];
        States(:,ctx+1) = x_step;
        % evaluate objective function
        TMus = sum(dMs(:,i).*FT);
        Terror = Tid(i) - TMus;
        J(ctx) = (Terror/10)^2;
        
        % aux variables
        TMusV(ctx) = TMus;
        
        % append counter states
        ctx = ctx+1;
    end
    if strcmp(Set.GainsEstimation,'GeyerDefault_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2') || ...
            strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
        f_Int1000_Swing_Full = Function('f_Int1000_Swing_Full',{OptReflex,ConstReflex,x0,lMT,dMs,Tid},{States,J,TMusV});
        f_Int1000_Swing_J = Function('f_Int1000_Swing_J',{OptReflex,ConstReflex,x0,lMT,dMs,Tid},{J});
    else
        f_Int1000_Swing_Full = Function('f_Int1000_Swing_Full',{OptReflex,x0,lMT,dMs,Tid},{States,J,TMusV});
        f_Int1000_Swing_J = Function('f_Int1000_Swing_J',{OptReflex,x0,lMT,dMs,Tid},{J});
    end
    f_Int1000_Swing_Full.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_Swing_Full']));
    f_Int1000_Swing_J.save(fullfile(OutPath,[Set.GainsEstimation 'f_Int1000_Swing_J']));    
    
    disp(['Create casadi functions finished for ' num2str(iSet) '/' num2str(length(Sets))]);
    clearvars -except iSet Sets MainPath OutPath
end

%% save casadifunctions



