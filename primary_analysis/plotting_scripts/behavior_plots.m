
%% BEHAVIOR CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Lucas Encarnacion-Rivera
%please consider citing the github <insert github here> if you found this
%useful 
%Description: this code plots a series of relevant behavioral metrics from extracted
%behavioral data during free foraging for food and water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load('\path\to\the\data')

%% reward collection over time 
figure;stairs(vertcat(rewcol.food,syncs(end)), [1:numel(rewcol.food)+1]); 
hold on;stairs(vertcat(rewcol.water,syncs(end)), [1:numel(rewcol.water)+1])
legend('food','water')
xlabel('time (s)');ylabel('rewards collected')

%% plot licking behavior (lick raster) 
foodcolor=[1 0.9 0.3];
watercolor=[0.2 0.7450 0.9330];
figure; hold on;y=[];
for x=1:sum(horzcat(length(order.order2),length(order.order4),length(order.order5),length(order.order9)))
    if portorder(1,x)==2 & length(order.order2) >= portorder(2,x)
        y(1:numel(order.order2{1,portorder(2,x)}))=x;
        s=scatter(order.order2{1,portorder(2,x)}/100,(y),10,foodcolor,'filled','o','AlphaData',.5);
        clear y;s.MarkerFaceAlpha=.75;
    end 
    if portorder(1,x)==4 & length(order.order4) >= portorder(2,x)
         y(1:numel(order.order4{1,portorder(2,x)}))=x;
        s=scatter(order.order4{1,portorder(2,x)}/100,(y),10,foodcolor,'filled','o','AlphaData',.5);
        clear y;s.MarkerFaceAlpha=.75;
    end 
    if portorder(1,x)==5 & length(order.order5) >= portorder(2,x)
        y(1:numel(order.order5{1,portorder(2,x)}))=x;
        s=scatter(order.order5{1,portorder(2,x)}/100,(y),10,watercolor,'filled','o','AlphaData',.5);
        clear y;s.MarkerFaceAlpha=.75;
    end 
    if portorder(1,x)==9 & length(order.order9) >= portorder(2,x)
         y(1:numel(order.order9{1,portorder(2,x)}))=x;
        s=scatter(order.order9{1,portorder(2,x)}/100,(y),10,watercolor,'filled','o','AlphaData',.5);
        clear y;s.MarkerFaceAlpha=.75;
    end 
end
hold off;
ylabel('reward bout');xlabel('time (s)');
xlim([0 1.5]);


%% All trajectories uniform orientation colored by consummatory outcome 

figure;hold on;
alphaval = 0.2;
for x = 1:numel(FZ.sames)
    nextsame = plot(FZ.sames{1,x}(:,1),FZ.sames{1,x}(:,2),'Color',[1, 0, 0, alphaval]);
end

xlabel('x-values (pixels)')
ylabel('y-values (pixels)')
title('Overlaid Trajectories')

alphaval = 0.5;
for x = 1:numel(FZ.switches)
    nextswitch = plot(FZ.switches{1,x}(:,1),FZ.switches{1,x}(:,2),'Color',[0, 0, 1, alphaval]);
end

alphaval = .8;
for x = 1:numel(FZ.misses)
    nextmiss = plot(FZ.misses{1,x}(:,1),FZ.misses{1,x}(:,2),'Color',[0, 1, 0, alphaval],'LineWidth',1);
end
plot(ploc(:,1),ploc(:,2),'d','LineWidth',3,'MarkerSize',5,'MarkerEdgeColor','k')
hold off

legend([nextsame,nextswitch,nextmiss], 'Same','Switch','Miss')

%% plot syllables in space 
figure;hold on;
for x=[2,3,4,5,6,7,1]
    scatter(fztrack.tracks{1,3}(1,find(kpmsfz.syllreindexed==x)),...
        fztrack.tracks{1,3}(2,find(kpmsfz.syllreindexed==x)), ...
        'SizeData',5,...
        'MarkerEdgeColor',kpmsfz.colors(x,:),'MarkerFaceColor',kpmsfz.colors(x,:),...
        'MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1);
end 

%% Plot ethogram of behavior over bouts

figure;hold on;xlim([10 300]);
for x=1:rewcol.boutcount
    g=repmat(x,numel(kpmsfz.syllbybout_cs{1,x}),1)';
    for i=kpmsfz.sylls'
        plot(find(kpmsfz.syllbybout_cs{1,x}==i),g(find(kpmsfz.syllbybout_cs{1,x}==i)),...
             'MarkerEdgeColor',kpmsfz.colors(i,:),'MarkerFaceColor',kpmsfz.colors(i,:),...
             'Marker','_','LineStyle','none','MarkerSize',1);
    end 
end 
xlabel('time (frames)');ylabel('reward sequence')
