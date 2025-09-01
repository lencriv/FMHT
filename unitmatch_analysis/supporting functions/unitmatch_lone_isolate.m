%% LOAD IN DATA
unitmatch_load

%% UNIQUE UNITS IN HUNGRY THIRSTY SESSIONS
%define um hungry thirsty vects 
uniqueIDs_ht=cell(1,length(mouseNames));
uniqueIDs_ht_sess=cell(1,length(mouseNames));
um_htvect{1}=1:5;um_htvect{2}=1:4;um_htvect{3}=1:5;um_htvect{4}=1:2;um_htvect{5}=1:3;
%loop over and append new neurons
for x = 1:length(mouseNames)
    uniqueIDs_ht{x}=[];
    for g=um_htvect{x} 
        h=unique(umdat{x}.MatchTable.UID1(umdat{x}.MatchTable.RecSes1==g));
        [~,i,~]=setxor(h,uniqueIDs_ht{x});
        uniqueIDs_ht{x}=vertcat(uniqueIDs_ht{x},h(i));
        uniqueIDs_ht_sess{x}{g}=h(i);
    end 
end 
%% computing thresholds
for x = 1:length(mouseNames)
    internalcheck.thresholdlone_mp(x)=prctile(umdat{x}.MatchTable.TotalScore(umdat{x}.MatchTable.ID1 == umdat{x}.MatchTable.ID2...
        & umdat{x}.MatchTable.RecSes1 == umdat{x}.MatchTable.RecSes2),3);
    internalcheck.thresholdlone_cd(x) = prctile(umdat{x}.MatchTable.CentroidDist(umdat{x}.MatchTable.ID1 == umdat{x}.MatchTable.ID2...
        & umdat{x}.MatchTable.RecSes1 == umdat{x}.MatchTable.RecSes2),50);
end 

%% identifying potential off matches 
for x=1:length(mouseNames)
    for g=1:length(uniqueIDs_ht_sess{x})
        h=0;
        for y=uniqueIDs_ht_sess{x}{g}' 
            h=h+1;
            internalcheck.matchprob{x}{g}{h} = umdat{x}.MatchTable.MatchProb(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);
            internalcheck.centroiddist{x}{g}{h} = umdat{x}.MatchTable.CentroidDist(umdat{x}.MatchTable.RecSes1==g &...
                        umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                            umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);
            internalcheck.UID2{x}{g}{h} = umdat{x}.MatchTable.UID2(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2<= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);
            internalcheck.recsess2{x}{g}{h} = umdat{x}.MatchTable.RecSes2(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2<= length(uniqueIDs_ht_sess{x}) &...
                         umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);
            internalcheck.matchprob_dist{x}{g}{h} = histcounts(internalcheck.matchprob{x}{g}{h},[0:.02:1],"Normalization","probability");
            internalcheck.matchprob_dist_all(h,:) = histcounts(internalcheck.matchprob{x}{g}{h},[0:.02:1],"Normalization","probability");
            internalcheck.matchprob_offmatchidx{x}{g}{h} = find(internalcheck.matchprob{x}{g}{h}>internalcheck.thresholdlone_mp(x) & ...
                internalcheck.matchprob{x}{g}{h}>internalcheck.thresholdlone_cd(x));
            internalcheck.UID_offmatch{x}{g}{h} = unique(internalcheck.UID2{x}{g}{h}(internalcheck.matchprob_offmatchidx{x}{g}{h}));
            internalcheck.UID_offmatch_rep{x}{g}{h} = internalcheck.UID2{x}{g}{h}(internalcheck.matchprob_offmatchidx{x}{g}{h});
            
            if ~isempty(internalcheck.UID_offmatch{x}{g}{h})
                internalcheck.offmatch_log{x}{g}(h) = 0;
            else 
                internalcheck.offmatch_log{x}{g}(h) = 1;
            end 
        end
        internalcheck.offmatch_count{x}(g)=sum(internalcheck.offmatch_log{x}{g});
    end 
end 

%% eliminating off matches for final set of unique units 
%get UIDs of correct matches using match log 
%for all the off matches, get the UID and the session and eliminate from
%the session with the greater neuron count 
k=0;
for x=1:length(mouseNames)
    for g=1:length(uniqueIDs_ht_sess{x})
        for y=find(internalcheck.offmatch_log{x}{g}==0)
            h=0;
            for k=internalcheck.UID_offmatch{x}{g}{y}'
                h=h+1;
                internalcheck.avg_offmatchprob{x}{g}{y}(h)=mean(internalcheck.matchprob{x}{g}{y}...
                    (internalcheck.matchprob_offmatchidx{x}{g}{y}(internalcheck.UID_offmatch_rep{x}{g}{y}==k)));
            end 
        end 
    end 
end 
%% algorithm to isolate lone neurons

for x=1:length(mouseNames)
    for g=1:length(uniqueIDs_ht_sess{x})
        h=0;d=0;
        internalcheck.final_lone_UIDs{x}{g}=[];
        internalcheck.disc_lone_UIDs{x}{g}=[];    
        for y=1:length(uniqueIDs_ht_sess{x}{g})
            if isempty(internalcheck.UID_offmatch{x}{g}{y}) 
                h=h+1;
                internalcheck.final_lone_UIDs{x}{g}(h)=uniqueIDs_ht_sess{x}{g}(y);
            elseif isscalar(internalcheck.UID_offmatch{x}{g}{y}) & ...
                   isempty(intersect(internalcheck.disc_lone_UIDs{x}{g},uniqueIDs_ht_sess{x}{g}(y)))
                d=d+1;h=h+1;
                internalcheck.final_lone_UIDs{x}{g}(h)=uniqueIDs_ht_sess{x}{g}(y);
                internalcheck.disc_lone_UIDs{x}{g}(d)=internalcheck.UID_offmatch{x}{g}{y};
            elseif length(internalcheck.UID_offmatch{x}{g}{y})>1 & isempty(intersect(internalcheck.final_lone_UIDs{x}{g},...
                    internalcheck.UID_offmatch{x}{g}{y}))
                %skip 
            end 
        end 
    end 
end 
%% get match probability for lone UIDs
%use final lone matches concatenated across sessions  
for x=1:length(mouseNames)
    for g=1:length(uniqueIDs_ht_sess{x})
        h=0;
        for y=internalcheck.final_lone_UIDs{x}{g} 
            h=h+1;
            internalcheck.matchprob_lone{x}{g}{h} = umdat{x}.MatchTable.MatchProb(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.tscore_lone{x}{g}{h} = umdat{x}.MatchTable.TotalScore(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.ISI_lone{x}{g}{h} = umdat{x}.MatchTable.ISICorr(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.centroiddist_lone{x}{g}{h} = umdat{x}.MatchTable.CentroidDist(umdat{x}.MatchTable.RecSes1==g &...
                        umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2 <= length(uniqueIDs_ht_sess{x}) &...
                            umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.UID2_lone{x}{g}{h} = umdat{x}.MatchTable.UID2(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2<= length(uniqueIDs_ht_sess{x}) &...
                        umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.recsess2_lone{x}{g}{h} = umdat{x}.MatchTable.RecSes2(umdat{x}.MatchTable.RecSes1==g &...
                    umdat{x}.MatchTable.RecSes2~=g & umdat{x}.MatchTable.RecSes2<= length(uniqueIDs_ht_sess{x}) &...
                         umdat{x}.MatchTable.UID1 == y & umdat{x}.MatchTable.UID2 ~= y);

            internalcheck.matchprob_dist_lone{x}{g}{h} = histcounts(internalcheck.matchprob_lone{x}{g}{h},[0:.02:1],"Normalization","probability");
            internalcheck.matchprob_dist_all_lone(h,:) = histcounts(internalcheck.matchprob_lone{x}{g}{h},[0:.02:1],"Normalization","probability");

            internalcheck.ISI_dist_lone{x}{g}{h} = histcounts(internalcheck.ISI_lone{x}{g}{h},[0:.02:1],"Normalization","probability");
            internalcheck.ISI_dist_all_lone(h,:) = histcounts(internalcheck.ISI_lone{x}{g}{h},[0:.02:1],"Normalization","probability");

            internalcheck.tscore_dist_lone{x}{g}{h} = histcounts(internalcheck.tscore_lone{x}{g}{h},[0:.02:1],"Normalization","probability");
            internalcheck.tscore_dist_all_lone(h,:) = histcounts(internalcheck.tscore_lone{x}{g}{h},[0:.02:1],"Normalization","probability");
            
            ID1_lone{x}{g}(h)=unique(umdat{x}.MatchTable.ID1(umdat{x}.MatchTable.RecSes1==g & umdat{x}.MatchTable.UID1==y));        
        end
    end 
end 

%% save variables
save(strcat(unitmatch_loc,'lone_neurons'),"internalcheck","ID1_lone");