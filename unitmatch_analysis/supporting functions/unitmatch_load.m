% clear;
npix_filepaths;
for nummice = 1:length(mouseNames)
    currentmouse = mouseNames{nummice};
    disp(strcat('loading',{' '},currentmouse,{' '},'match data'))
    load(strcat(unitmatch_loc,currentmouse,'\UnitMatch.mat'));
    MatchTable_old = MatchTable;
    load(strcat(unitmatch_loc,currentmouse,'\MatchingScores.mat'));
    load(strcat(unitmatch_loc,currentmouse,'\AUC.mat'));
    load(strcat(unitmatch_loc,currentmouse,'\UnitMatchModel.mat')); 
    load(strcat(unitmatch_loc,currentmouse,'\outputeval.mat')); 
    load(strcat(unitmatch_loc,currentmouse,'\clusinfo.mat')); 
    load(strcat(unitmatch_loc,currentmouse,'\',currentmouse,'_match_data.mat')); 
    load(strcat(unitmatch_loc,currentmouse,'\corrected_MatchTable.mat')); 

    session.clusinfo = clusinfo;session.MatchTable = MatchTable;session.clusinfo = clusinfo;
    session.AUCStruct = AUCStruct;session.BestMdl = BestMdl;session.clusinfo = clusinfo;session.WaveformInfo = WaveformInfo;
    session.GoodRecSesID = GoodRecSesID;session.matches_idx = matches_idx;session.matchesidx_all = matchesidx_all;
    session.MatchProbability = MatchProbability;session.MatchTable = MatchTable;session.modeleval = modeleval;
    session.neighbors_idx = neighbors_idx;session.samesessidx = samesessidx;session.self_idx = self_idx;
    session.TotalScore = TotalScore;session.UMparam = UMparam;session.uniqueID_recurrence = uniqueID_recurrence;
    session.uniqueID_session_presence = uniqueID_session_presence;session.UniqueIDConversion = UniqueIDConversion;
    session.uniqueIDs_all=uniqueIDs_all;session.uniqueID_recurrence=uniqueID_recurrence;session.threshold_lone=threshold_lone;
    session.uniqueID_session_presence=uniqueID_session_presence;session.MatchTable_old = MatchTable_old;
    umdat{nummice}=session;
    clearvars -except umdat htrecall trecall hrecall
    npix_filepaths
end 

%load session library
if exist(strcat(unitmatch_loc,'session_names_lib.mat'),"file") == 2
    load(strcat(unitmatch_loc,'session_names_lib'));
    for x=1:length(sessname_lib)
        for g=1:length(sessname_lib{x})
            sesstype_lib{x}(g)=session_type{find(vertcat(npixfilenames{:})==sessname_lib{x}{g})};
        end 
    end 
end 

%load lone neuron data
if exist(strcat(unitmatch_loc,'lone_neurons.mat'),"file") == 2
    load(strcat(unitmatch_loc,'lone_neurons'));
end 