
%% run unitmatch on all data 
clear
npix_filepaths
for um_counter=1:length(mouseNames)
    % insert a list of your filepaths for ephys output
    % insert a list of your filepaths for unitmatch
    UMparam.SaveDir = strcat('C:\Users\LEncR\Dropbox\motiv\FMHT_data\unitmatch\export\', ... 
                             mouseNames{um_counter});
    UMparam.KSDir = KSDir_um{um_counter};
    UMparam.RawDataPaths = RawDataPaths_um{um_counter};
    optimizeum=0;
    unitmatch_matching_code;  %this is modified from 'https://github.com/EnnyvanBeest/UnitMatch'
    clear;close all; 
    npix_filepaths;
end 

%% evaluate model 
unitmatch_evaluate_models

%% re-run with optimized parameters
clear
npix_filepaths
for um_counter=1:length(mouseNames)
    npix_filepaths;
    unitmatch_filepaths;
    UMparam.SaveDir = strcat('C:\Users\LEncR\Dropbox\motiv\FMHT_data\unitmatch\export\', ... 
                             mouseNames{um_counter});
    UMparam.KSDir = KSDir_um{um_counter}; 
    UMparam.RawDataPaths = RawDataPaths_um{um_counter};
    optimizeum=1; 
    unitmatch_matching_code; %this is modified from 'https://github.com/EnnyvanBeest/UnitMatch'
    clear;close all; 
    npix_filepaths;
end 

%% post processing 
unitmatch_postprocess
clear
%% isolate lone matches 
unitmatch_load 
unitmatch_lone_isolate




