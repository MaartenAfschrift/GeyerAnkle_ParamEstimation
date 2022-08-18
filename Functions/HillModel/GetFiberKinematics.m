function [cos_alpha, lM, lMT, vMT, dM, FL, FV, FP] = GetFiberKinematics(ModelPath,muscleNames,PolyPath,OptInfo)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

load(fullfile(PolyPath,'muscle_spanning_joint_INFO_CE.mat'));
load(fullfile(PolyPath,'MuscleInfo_CE.mat'));
load('ActiveFVParameters.mat','ActiveFVParameters');
load('PassiveFLParameters','PassiveFLParameters');
load('Faparam.mat','Faparam');

% read the muscle parameters
[params,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(ModelPath,muscleNames);

% get function to compute muscle geometry
NMuscle 	= length(MuscleInfo.muscle);
ndof        = 3;
qin     	= SX.sym('qin',1,ndof);
qdotin  	= SX.sym('qdotin',1,ndof);
lMT     	= SX(NMuscle,1);
vMT     	= SX(NMuscle,1);
dM      	= SX(NMuscle,ndof);
for i=1:NMuscle     
    index_dof_crossing  = find(muscle_spanning_joint_INFO(i,:)==1);
    order               = MuscleInfo.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
    lMT(i,1)            = mat*MuscleInfo.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:ndof)        = 0;
    nr_dof_crossing     = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*qdotin(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

% compute muscle geometry info
nfr = length(OptInfo.q);   cos_alpha = zeros(nfr,NMuscle); lM = zeros(nfr,NMuscle);
lMT = zeros(nfr,NMuscle);  vMT       = zeros(nfr,NMuscle); dM = zeros(nfr,NMuscle,ndof);
FL 	= zeros(nfr,NMuscle);  FV        = zeros(nfr,NMuscle); FP = zeros(nfr,NMuscle);
for ifr = 1:nfr
    [lMTs,vMTs,dMs] = f_lMT_vMT_dM(OptInfo.q(ifr,:),OptInfo.qd(ifr,:));
    lMT(ifr,:) = full(lMTs);    vMT(ifr,:) = full(vMTs);    dM(ifr,:,:) = full(dMs);       
end

for i=1:length(lMT(:,1))
    [FL(i,:),FV(i,:),FP(i,:),cos_alpha(i,:),lM(i,:)] = ...
        HillModel_RigidTendon(lMT(i,:),vMT(i,:),params,ActiveFVParameters,PassiveFLParameters,Faparam);
end



end

