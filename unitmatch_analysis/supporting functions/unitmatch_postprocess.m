clear;
npix_filepaths;
for nummice = 1:length(mouseNames)
    %%%%% load unitmatch data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    npix_filepaths;
    currentmouse = mouseNames{nummice};
    
    load(strcat(unitmatch_loc,currentmouse,'\UnitMatch.mat'));
    load(strcat(unitmatch_loc,currentmouse,'\MatchingScores.mat'));
    load(strcat(unitmatch_loc,currentmouse,'\AUC.mat'));
    load(strcat(unitmatch_loc,currentmouse,'\UnitMatchModel.mat'));

    disp(strcat('matching',{' '},currentmouse))
    %%%%% initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    samesessidx=cell(1, length(mousesessidx{nummice}));
    neighbors_idx=cell(1, length(mousesessidx{nummice}));
    self_idx=cell(1, length(mousesessidx{nummice}));  
    matches_idx=cell(1, length(mousesessidx{nummice}));
    goodUniqueIDs = unique(UniqueIDConversion.UniqueID(find(UniqueIDConversion.GoodID)));
    uniqueID_recurrence = zeros(1,length(goodUniqueIDs));
    uniqueID_session_presence= zeros(length(mousesessidx{nummice}),length(goodUniqueIDs));
    %%%%% match each session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load in data
    for mousesessnum = 1:length(mousesessidx{nummice})
        um.currsessionname = npixfilenames{1,mousesessidx{nummice}(mousesessnum)};
        
        disp(strcat('loading',{' '},um.currsessionname,{' '},num2str(mousesessnum),{' '},'of',{' '},num2str(length(mousesessidx{nummice}))))        
        %map original cluster ID from session onto unitmatch cluster IDs 
        
        %sessidx in unitmatch 
        for x=unique(GoodRecSesID)'
            if contains(UniqueIDConversion.Path4UnitNPY(find(GoodRecSesID==x,1)),um.currsessionname) == 1 
               um.currsessumidx = x;
            end     
        end 
        % get unique IDs of units in session 
        sessname_lib{nummice}{um.currsessumidx}=um.currsessionname; 
        um.sessidx = find(GoodRecSesID == um.currsessumidx); 
        um.uniqueIDs = unique(MatchTable.UID1(MatchTable.RecSes1==um.currsessumidx));
        um.clusIDs = UniqueIDConversion.OriginalClusID(um.sessidx);

        %% fix redundant matches and reassign UID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        unitmatch_fix_redundancies
        
        %% new IDs 
        um.uniqueIDs = unique(MatchTable.UID1(MatchTable.RecSes1==um.currsessumidx));
        um.clusIDs = unique(MatchTable.ID1(MatchTable.RecSes1==um.currsessumidx));

        %%%%% get neighbor and match distributions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        samesessidx{mousesessnum} = find(MatchTable.RecSes1-MatchTable.RecSes2 == 0 & ...
                                         MatchTable.RecSes1+MatchTable.RecSes2 == um.currsessumidx * 2);
                                         %^index of same session comparisons
        self_idx{mousesessnum} = intersect(samesessidx{mousesessnum},find(MatchTable.ID1 - MatchTable.ID1 == 0));
                    %^from the same session get index where the same units are being compared (same distribution)
        matches_idx{mousesessnum} = find(MatchTable.RecSes1==um.currsessumidx & ...
                                         MatchTable.UID1 == MatchTable.UID2 & ... 
                                         MatchTable.RecSes1 ~= MatchTable.RecSes2); 
                    %^from the same session get index where matched units are being compared (match distribution)
        neighbors_idx{mousesessnum} = intersect(samesessidx{mousesessnum},...
                                      find(MatchTable.ID1 ~= MatchTable.ID2));
                    %^from the same session get index where the different units are being compared (neighbor distribution)
    

        um.samesessidx=samesessidx{mousesessnum};
        um.self_idx=self_idx{mousesessnum};
        um.neighbors_idx=neighbors_idx{mousesessnum};
        um.matches_idx=matches_idx{mousesessnum};
        %% get match probabilities 
        um.matchUIDs = unique(MatchTable.UID1(um.matches_idx));
        um.matchclusIDs = unique(MatchTable.ID1(um.matches_idx));
        h=0;um.match_matchprobs=[];
        for x=um.matchUIDs'
            h=h+1; 
            um.match_matchprobs{h} = MatchTable.MatchProb...
                (find(MatchTable.UID1 == x & MatchTable.RecSes1 == um.currsessumidx)); 
        end 

        %% get lone IDs 
        um.loneUIDs = setxor(um.uniqueIDs,um.matchUIDs); 
        h=0;
        for x=um.loneUIDs' 
            h=h+1;
            um.lone_matchprobs{h} = MatchTable.MatchProb(MatchTable.UID1==x & ...
                MatchTable.RecSes1 ~= MatchTable.RecSes2); 
            um.loneclusIDs(h) = unique(MatchTable.ID1(MatchTable.UID1==x));
        end 

        threshold_lone = .85 - 2 * std(MatchTable.MatchProb(um.matches_idx));

        um.lone_log = zeros(1,length(um.lone_matchprobs));
        for x=1:length(um.lone_matchprobs)
            if isempty(find(um.lone_matchprobs{x}>threshold_lone))
                um.lone_log(x)=1;
            end 
        end 
        %% SAVE OUT UM STRUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(strcat('saving',{' '},um.currsessionname,{' '},num2str(mousesessnum),{' '},'of',{' '},num2str(length(mousesessidx{nummice}))))
        save(strcat('Z:\motiv\data\export\',um.currsessionname,'_um'), 'um')
        clear um isredundant
    end 
    matchesidx_all = find(table2array(MatchTable(:,5)) - table2array(MatchTable(:,6)) == 0 & ... 
                      table2array(MatchTable(:,3)) - table2array(MatchTable(:,4)) ~= 0); 
                        % ^ index of comparisons of matches between sessions
    %% get recurrence of matched IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    uniqueIDs_all=unique(MatchTable.UID1);
    uniqueID_session_presence = zeros(length(mousesessidx{nummice}),length(uniqueIDs_all));

    for g=1:length(mousesessidx{nummice})
        [~,~,i]=intersect(unique(MatchTable.UID1(MatchTable.RecSes1==g)),uniqueIDs_all);
        uniqueID_session_presence(g,i')=1;
    end 

    uniqueID_recurrence = sum(uniqueID_session_presence);
    %%%%% SAVE OUT ALL MATCHTABLE IDXs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(strcat('saving',{' '},currentmouse,{' '},num2str(mousesessnum),{' '},'of',{' '},num2str(length(mousesessidx{nummice}))))
    %saveout corrected match table 
    save(strcat(unitmatch_loc,'\',currentmouse,'\corrected_MatchTable'),"MatchTable");
    %saveout match data 
    save(strcat(UMparam.SaveDir,'\',currentmouse,'_match_data'), "self_idx","neighbors_idx","uniqueIDs_all",...
                                                 "matchesidx_all","matches_idx","samesessidx","MatchTable",...
                                                 "uniqueID_recurrence","uniqueID_session_presence","threshold_lone")
    %saveout session idx
    save(strcat(UMparam.SaveDir,'\' ,'session_names_lib'),"sessname_lib"); 
    clear
end 
