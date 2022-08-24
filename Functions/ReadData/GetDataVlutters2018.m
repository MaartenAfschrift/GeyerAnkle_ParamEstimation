function [DataPath] = GetDataVlutters2018()
%GetDataVlutters2018 Downloads the data of the pelvis push and pull
%perturbations. This dataset was organized in such a way that it is
%convenient for the parameter estimation.

% Publication Vlutters 2018: 
% Vlutters, M., van Asseldonk, E.H.F. & van der Kooij, H. 
% Lower extremity joint-level responses to pelvis perturbation during human
% walking. Sci Rep 8, 14621 (2018). https://doi.org/10.1038/s41598-018-32839-8
% 
% get path to DataFolder
[ReadDataPath,~,~] = fileparts(mfilename('fullpath'));
[FuncPath,~] = fileparts(ReadDataPath);
[MainPath,~] = fileparts(FuncPath);
DataPath = fullfile(MainPath,'Data');

if isfolder(DataPath) && isfile(fullfile(DataPath,'pp_1','PertWalk.mat'))
    disp('Data was already downloaded');
else
    if ~isfolder(DataPath)
        mkdir(DataPath);
    end
    % download the executable from google drive
    disp('Start Downloading vlutters data (320 MB)');
    LinkZipFile = 'https://www.dropbox.com/s/c0vl9skfo08v7ig/Vlutters2018_ParamID.zip?dl=1';    
    websave(fullfile(DataPath,'Vlutters2018_ParamID.zip'),LinkZipFile);
    disp('unzip DataFile');
    unzip(fullfile(DataPath,'Vlutters2018_ParamID.zip'),DataPath);
    delete(fullfile(DataPath,'Vlutters2018_ParamID.zip'));
    disp('Downloading Vlutters 2018 data finished');
end

end