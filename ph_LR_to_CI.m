function pop=ph_LR_to_CI(keys,pop)
%% assigning contra and ipsi instead of left/right, dependent on the recorded target (which contains hemisphere information in the end)
% reference is left hemisphere, meaning that right becomes contra and left ipsi
% that is why hemifields, positions and hands need to be inverted for right hemisphere targets.
% for right recording sites LR->CI  (hemifield & positions...?)
% positive space and hand==2 become contra
% see https://github.com/dagdpz/spike_analysis/blob/master/ph_LR_to_CI.m
% borrowed from spike_analysis

target=pop.(keys.contra_ipsi_relative_to);
if strcmpi(target(end-1:end),'_R') || strcmpi(target(end-1:end),'_r')
    poptr=num2cell([pop.trial.hemifield]*-1);
    [pop.trial.hemifield]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{pop.trial.position},'uniformoutput',false);
    [pop.trial.position]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{pop.trial.fixation},'uniformoutput',false);
    [pop.trial.fixation]=deal(poptr{:});
    hand1=[pop.trial.reach_hand]==2;
    hand2=[pop.trial.reach_hand]==1;
    [pop.trial(hand1).reach_hand]=deal(1);
    [pop.trial(hand2).reach_hand]=deal(2);
end
end
