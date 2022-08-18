function [OptReflex, x0] = UnpackSolutionShooting_AnkleOnly(w_opt,Set,np)
% UnpackSolutionShooting_AnkleOnly Function to unpack the w_opt vector after a shooting simulation
% with the geyer model with only the soleus and tibialis anterior muscles

% get the initial state
if Set.OptInitState
	x0 = zeros(6,np);
	for i=1:np
	    x0(:,i) = w_opt(i*6-5:i*6);
	end

	% number of optimization variables for the reflexes depend on the case used
	i0 = i*6+1;
else
    i0 = 1;
    x0 = NaN;
end

% get the reflex parameters
if strcmp(Set.GainsEstimation,'GeyerDefault') || strcmp(Set.GainsEstimation,'Geyer_Default_GRF')
    OptReflex = w_opt(i0:i0+6-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COM_COMd')
    OptReflex = w_opt(i0:i0+10-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd') || strcmp(Set.GainsEstimation,'Geyer_COMd_GRF')
    OptReflex = w_opt(i0:i0+8-1);
elseif strcmp(Set.GainsEstimation,'GeyerDefault_Limited1')  
    OptReflex = w_opt(i0:i0+3-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited2')
    OptReflex = w_opt(i0:i0+5-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_Limited1')
    OptReflex = w_opt(i0:i0+6-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMd_GRF_lim')
    OptReflex = w_opt(i0:i0+7-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMAll_GRF')
    OptReflex = w_opt(i0:i0+12-1);
elseif strcmp(Set.GainsEstimation,'Geyer_COMdd_GRF')
    OptReflex = w_opt(i0:i0+10-1);
end

end


