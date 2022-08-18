function [Data] = GetMuscleModel_vTwente(PolyPath,ModelPath,muscleNames,Data,headers,varargin)
%GetMuscleModel_vTwente Get all the hilltype muscle parameters and geometry
%   Detailed explanation goes here

BoolConvertRad = 0;
if ~isempty(varargin)
    BoolConvertRad = varargin{1};
end

% offset in angle angle
qa_offset = 0;
if length(varargin)>1
    qa_offset = varargin{2};
end

import casadi.*
load(fullfile(PolyPath,'muscle_spanning_joint_INFO_CE.mat'));
load(fullfile(PolyPath,'MuscleInfo_CE.mat'));

% read the muscle parameters
[params]=ReadMuscleParameters(ModelPath,muscleNames);

% loop over all trials
ntr = length(Data);
for itr = 1:ntr
%     disp(num2str(itr))
    % evaluate spline
    tspan = Data(itr).tspan;
    tSample = tspan(1):0.01:tspan(end);
    DatSpline = ppval(Data(itr).DataSpline,tSample)';
    
    % evaluate derivative
    SplineSel_dot = fnder(Data(itr).DataSpline,1);
    DatSpline_dot = ppval(SplineSel_dot,tSample)';
    
    iqa= find(strcmp(headers,'qa'));
    iqk= find(strcmp(headers,'qk'));
    iqh= find(strcmp(headers,'qh'));
    q = DatSpline(:,[iqa iqk iqh]);
    qd = DatSpline_dot(:,[iqa iqk iqh]);
    
    NMuscle = length(MuscleInfo.muscle);
    ndof = 3;
    
    % convert to radians if needed
    if BoolConvertRad
        q = q * pi / 180;
        qd = qd * pi / 180;
    end
    
    % offset in angles
    q(:,1) = q(:,1) + qa_offset;
    
    
    % adapt q's to definition polynomials - >[hip knee ankle]    
    q = [q(:,3) q(:,2) q(:,1)];
    qd = [qd(:,3) qd(:,2) qd(:,1)];
    
    
    % compute the muscle-tendon length and velocity based on the polynomials
    nfr = length(tSample); lMT = zeros(nfr,NMuscle);
    vMT = zeros(nfr,NMuscle); dM = zeros(nfr,NMuscle,ndof);
    for i=1:NMuscle
        index_dof_crossing  = find(muscle_spanning_joint_INFO(i,:)==1);
        order               = MuscleInfo.muscle(i).order;
        [mat,diff_mat_q]    = n_art_mat_3(q(:,index_dof_crossing),order);
        nr_dof_crossing     = length(index_dof_crossing);
        lMT(:,i)            = mat*MuscleInfo.muscle(i).coeff;
        for dof_nr = 1:nr_dof_crossing
            dM(:,i,index_dof_crossing(dof_nr)) = (-(diff_mat_q(:,:,dof_nr)))*MuscleInfo.muscle(i).coeff;
            vMT(:,i) = vMT(:,i) + (-dM(:,i,index_dof_crossing(dof_nr)).*qd(:,index_dof_crossing(dof_nr)));
        end
    end
    
    % spline fit of muscle-tendon length and velocity
    Data(itr).lMT_Spline = spline(tSample,lMT');
    for dof = 1:ndof
        Data(itr).dof(dof).dM_Spline = spline(tSample,squeeze(dM(:,:,dof))');
    end
    Data(itr).vMT_Spline = spline(tSample,vMT');
    
    % Store Muscle Info
    Data(itr).Muscle.params = params;
    Data(itr).Muscle.Nmuscle = NMuscle;
end


end

