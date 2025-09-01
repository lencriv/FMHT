function plot_centered_activity_bout_fw(neuronidx,FZ, behave, spikes_neuron,recording_info,consumeendd,consume,region,htrecsess)
    fidx=sort(intersect(find(FZ.allrews(2,:)==2 | FZ.allrews(2,:)==4),...
        (horzcat(intersect(behave.reorient.left_boutidx,behave.approach.boutidx),...
                 intersect(behave.reorient.right_boutidx,behave.approach.boutidx)))));
    widx=sort(intersect(find(FZ.allrews(2,:)==5 | FZ.allrews(2,:)==9),...
        (horzcat(intersect(behave.reorient.left_boutidx,behave.approach.boutidx),...
                 intersect(behave.reorient.right_boutidx,behave.approach.boutidx)))));
    foodcolor=[1 0.9 0.3];
    watercolor=[0.2 0.7450 0.9330];
    figure;subplot(3,1,2);hold on;;h=0;k=0;j=[];
    for g=randfromvect(widx,100)
        k=k+1;
        l=(spikes_neuron{neuronidx}(spikes_neuron{neuronidx}>recording_info.frames_binned_fz(consumeendd(g)) - .5 & ...
            spikes_neuron{neuronidx}<recording_info.frames_binned_fz(consumeendd(g)) + 3)) -...
                recording_info.frames_binned_fz(consumeendd(g));
        j=vertcat(j,histcounts(l,[-0.5:0.05:2]));
        scatter(l,repmat(k,1,numel(l)),1,'MarkerEdgeColor',watercolor,'MarkerFaceColor','k');
        scatter(recording_info.frames_binned_fz(consume(1,g+1))-...
            recording_info.frames_binned_fz(consume(2,g)),k,20,'r','filled','square')
        xlim([-0.5,3])
    end 
    xlim([-0.5,3]);xline(0,'r');subplot(3,1,1);hold on
    plot([-0.5:0.05:1.95],mean(j),'Color',watercolor);xlim([-0.5,2]);xline(0,'r');
    %food
    subplot(3,1,3);hold on;h=0;k=0;j=[];
    for g=randfromvect(fidx,100)
        k=k+1;
        l=(spikes_neuron{neuronidx}(spikes_neuron{neuronidx}>recording_info.frames_binned_fz(consumeendd(g)) - .5 & ...
            spikes_neuron{neuronidx}<recording_info.frames_binned_fz(consumeendd(g)) + 3)) -...
                recording_info.frames_binned_fz(consumeendd(g));
        j=vertcat(j,histcounts(l,[-0.5:0.05:2]));
        scatter(l,repmat(k,1,numel(l)),1,'MarkerEdgeColor',foodcolor,'MarkerFaceColor','k');
        scatter(recording_info.frames_binned_fz(consume(1,g+1))-...
            recording_info.frames_binned_fz(consume(2,g)),k,20,'r','filled','square')
        xlim([-0.5,3])
    end 
    xlim([-0.5,3]);xline(0,'r');subplot(3,1,1);hold on
    plot([-0.5:0.05:1.95],mean(j),'Color',foodcolor);xlim([-0.5,2]);xline(0,'r');
    title(strcat(table2array(region.names(neuronidx,1)),'neuron',num2str(neuronidx),'sessid',num2str(htrecsess)));
end 