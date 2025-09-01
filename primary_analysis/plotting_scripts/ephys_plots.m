
%% NEURAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Lucas Encarnacion-Rivera
%please consider citing the github <insert github here> if you found this
%useful 
%Description: this code plots a series of relevant plots on neural data
% from neuropixels recording along a trajectory spanning hippocampus,
% thalamus, and hypothalamus, while an animal is freely foraging for food
% and water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Whole session firing rates with reward and room location

%plots the neural activity of a whole representative session with the times
%of food and water collection and the zone that the mouse was in 

g=find(~isnan(fztrack.tracks{6}(1,:)));h=find(isnan(fztrack.tracks{6}(1,:)));
figure;subplot(3,1,1);hold on;xlabel('time');ylabel(' ');
l=gca; l.XAxis.TickLength = [0 0]; 
set(gca, 'Position', [0.1 0.82 0.8 0.12]); % Adjusted to be thin

for x=1:50;scatter(h,repmat(x,1,length(h)),10,'red','filled','square');end 
for x=1:50;scatter(g,repmat(x,1,length(g)),10,'green','filled','square');end 
hold off;xlim([0,trkdim(1)]);

subplot(3,1,2);hold on;l=gca; l.XAxis.TickLength = [0 0]; 
set(gca, 'Position', [0.1 0.62 0.8 0.12]); % Adjusted to be thin

g=consume(1,FZ.allrews(2,:)==2 | FZ.allrews(2,:)==4);
h=consume(1,FZ.allrews(2,:)==5 | FZ.allrews(2,:)==9);

for x=1:50;scatter(g,repmat(x,1,length(g)),10,foodcolor,"filled","square");end  
for x=1:50;scatter(h,repmat(x,1,length(h)),10,watercolor,"filled","square");end 
xlim([0,trkdim(1)]);

subplot(3,1,3);set(gca, 'Position', [0.1 0.07 0.8 0.5]);
imagesc(zfr_matrix_clustered);
caxis([-1 3]);colormap('parula');c = colorbar('Location', 'southoutside');ylabel('neurons');

%% Example neurons 

%plots the firing rates with raster plots of single neuron activity  
%centered on bout start, and different actions such as reorient left,
%right, and approach.

%single neuron activity centered on reorient right
plot_centered_activity_rr(61,recording_info,behave,spikes_neuron,region,6)
%single neuron activity centered on reorient left  
plot_centered_activity_rl(89,recording_info,behave,spikes_neuron,region,6)
%single neuron activity centered on approach
plot_centered_activity_approach(40,recording_info,behave,spikes_neuron,region,6)
%single neuron activity centered on bout start separated by food and water
plot_centered_activity_bout_fw(66,FZ, behave, spikes_neuron,recording_info,consumeendd,consume,region,6)

%% projection onto drive dimension

figure;subplot(2,1,1)
linearprojlda(ldatest.same.ldaModel,recording_info.noconsume_frps(:,fwidxs.fwidx),sampleact.same.fwidx);
title('food vs water');
xlim([0 length(fwidxs.fwidx)]);ylim([-0.42 -0.3])



