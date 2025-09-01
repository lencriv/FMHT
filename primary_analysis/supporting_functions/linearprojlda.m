
%% LINEAR PROJECTION OF LDA COEFFICIENTS
function linearprojlda(model,feats_by_trials,trialidx) 
    coefficients = model.Coeffs(1, 2).Linear; % For two classes, adjust as needed for more classes
    trials = 1:size(feats_by_trials, 2);
    bias = model.Coeffs(1, 2).Const
    lda_projection = (feats_by_trials' * coefficients)+bias; % Linear projection (trials x 1)
    lda_projection = negonetoonenorm(lda_projection);
    plot(trials, smoothdata(lda_projection,'movmean',5),'Color','k','LineStyle', '-','LineWidth',2); 
    xlabel('Trial');
    ylabel('LDA coefficients (arbitrary units)');
    title('LDA Projection of Trials');
    hold on;
    scatter(trials(trialidx == 1), lda_projection(trialidx == 1), 'r', 'filled','MarkerFaceAlpha',.5); % Class 1 in red
    scatter(trials(trialidx == 0), lda_projection(trialidx == 0), 'b', 'filled','MarkerFaceAlpha',.5); % Class 0 in blue
    legend('smoothed mean', 'Class 1', 'Class 0');
    hold off;
end 