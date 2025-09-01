%% Primary pipeline to generate figures from example neural and behavioral data 
%Author: Lucas Encarnacion-Rivera 
%this script generates plots related to neural and behavioral data from a
%mouse foraging for food and water. 

%% Load in the data      
load('\path\to\the\data')
%% add path to access functions
addpath('\path\to\the\functions');
%% behavior plots
behavior_plots
%% neural data plots 
ephys_plots 
