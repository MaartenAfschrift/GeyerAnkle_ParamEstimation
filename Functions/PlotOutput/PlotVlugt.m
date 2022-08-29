function [h] = PlotVlugt(P,varargin)
%PlotVlugt Visualises the co-variance matrix of the jacobian of the cost
%function w.r.t. the parameters
%   P: covariance matrix
%   varargin:
%       (1) string with ID of the parameter set (needed for labels)

[nr,nc] = size(P);

if isempty(varargin)
    if nr == 8
        % assumed full set of gains
        Labels ={'Sol e0','Sol G','Tib e0','Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};
    elseif nr == 6
        % assumed default geyer model
        Labels ={'Sol e0','Sol G','Tib e0','Tib G','Tib loff','Tib Gsol'};
    elseif nr == 3
        % assumed limited set default
        Labels ={'Sol G','Tib G','Tib Gsol'};
    elseif nr == 5
        % assumed limited set with COM fb
        Labels ={'Sol G','Tib G','Tib Gsol','Sol COMd','Tib COMd'};
    elseif nr == 7
        % assumed limited set with COM fb
        Labels ={'Sol G','Tib e0','Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};
    else
        error('Unkown set of feedback gains');
    end
    ParamID = [];
else
    ParamID = varargin{1};
    if strcmp(ParamID,'GeyerDefault')
        Labels ={'Sol e0','Sol G','Tib e0','Tib G','Tib loff','Tib Gsol'};        
    elseif strcmp(ParamID,'Geyer_COM_COMd')
        Labels ={'Sol e0','Sol G','Tib e0','Tib G','Tib loff','Tib Gsol','Sol COM','Sol COMd','Tib COM','Tib COMd'};
    elseif strcmp(ParamID,'Geyer_COMd') || strcmp(ParamID,'Geyer_COMd_GRF')
        Labels ={'Sol e0','Sol G','Tib e0','Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};
    elseif strcmp(ParamID,'GeyerDefault_Limited1')
        Labels ={'Sol G','Tib G','Tib Gsol'};
    elseif strcmp(ParamID,'Geyer_COMd_Limited2')
        Labels ={'Sol G','Tib G','Tib Gsol','Sol COMd','Tib COMd'};
    elseif strcmp(ParamID,'Geyer_COMd_Limited1')
        Labels ={'Sol G','Tib e0','Tib G','Tib Gsol','Sol COMd','Tib COMd'};
    elseif strcmp(ParamID,'Geyer_COMd_GRF_lim')
        Labels ={'Sol G','Tib e0','Tib G','Tib loff','Tib Gsol','Sol COMd','Tib COMd'};
    end
    
end
nLabels = length(Labels);


h = figure('Name',ParamID);
image(abs(P),'CDataMapping','scaled');
set(gca,'CLim',[0 1]);
colorbar;

set(gca,'XTick',1:nLabels);
set(gca,'YTick',1:nLabels);
set(gca,'XTickLabel',Labels);
set(gca,'YTickLabel',Labels);

title(['Interdependence parameters: ' ParamID]);


delete_box;
end

