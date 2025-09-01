clear
npix_filepaths
close all
for nummice = 1:length(mouseNames)
    currentmouse = mouseNames{nummice}
    EvaluatingUnitMatch(strcat('C:\Users\LEncR\Dropbox\motiv\FMHT_data\unitmatch\export\',currentmouse)); 
    figure(1); 
    saveas(gcf, fullfile(strcat('C:\Users\LEncR\Dropbox\motiv\FMHT_data\unitmatch\export\',...
        currentmouse,'_units_tracked'))); 
    figure(2); 
    saveas(gcf, fullfile(strcat('C:\Users\LEncR\Dropbox\motiv\FMHT_data\unitmatch\export\',...
        currentmouse,'\',currentmouse,'_unit_parameters_distributions'))); 
    close all
end 
