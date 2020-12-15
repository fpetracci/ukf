%This script loads in the Matlab path the folder cointaining the ukf function.

%%
disp('UKF: Alessia Biondi and Francesco Petracci')

if verLessThan('matlab', '7.0')
    warning('You are running a very old (and unsupported) version of MATLAB.  You will very likely encounter significant problems using this files but you are on your own with this');
end


script_path = fileparts( mfilename('fullpath') );
addpath(script_path);

%addpath(genpath(fullfile(script_path, 'functions')));

clear script_path