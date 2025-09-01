function plot_centered_activity_rr(neuron_idx,recording_info,behave,spikes_neuron,region,htrecsess)
    k=0;j=[];
    figure;subplot(2,1,2);hold on;ylabel('bout');xlabel('time (s)')
    for g=randfromvect([1:length(behave.reorient.right_starts)],100)
        k=k+1;
        l=(spikes_neuron{neuron_idx}(spikes_neuron{neuron_idx}>recording_info.frames_binned_fz(behave.reorient.right_starts(g)) - .5 & ...
            spikes_neuron{neuron_idx}<recording_info.frames_binned_fz(behave.reorient.right_starts(g)) + 2)) -...
                recording_info.frames_binned_fz(behave.reorient.right_starts(g));
        j=vertcat(j,histcounts(l,[-0.5:0.05:2]));
        scatter(l,repmat(k,1,numel(l)),1,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(recording_info.frames_binned_fz(behave.reorient.right_ends(g))-...
            recording_info.frames_binned_fz(behave.reorient.right_starts(g)),k,20,'r','filled','square')
    end 
    xlim([-0.5,2]);xline(0,'r');
    subplot(2,1,1);plot([-0.5:0.05:1.95],mean(j),'k');xlim([-0.5,2]);xline(0,'r');
    title(strcat(table2array(region.names(neuron_idx,1)),'neuron',num2str(neuron_idx),'sessid',num2str(htrecsess)));
    ylabel('density')
end 