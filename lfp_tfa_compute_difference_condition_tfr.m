function [ diff_tfr ] = lfp_tfa_compute_difference_condition_tfr( lfp_tfr, diff_condition, stat_test, lfp_tfa_cfg )
%lfp_tfa_compute_difference_condition_tfr - function to compute the difference
%of LFP time frequency spectrogram averages between different conditions
%
% USAGE:
%   diff_tfr = lfp_tfa_compute_difference_condition_tfr( lfp_tfr,
%	diff_condition)
%	diff_tfr = lfp_tfa_compute_difference_condition_tfr( lfp_tfr,
%	diff_condition, stat_test, lfp_tfa_cfg )
%
% INPUTS:
%       lfp_tfr   - struct containing the LFP time frequency power
%       spectrogram averages for different conditions as average across
%       multiple trials within a site, or average across sites within
%       single or multiple sessions, see lfp_tfa_plot_site_average_tfr,
%       lfp_tfa_avg_tfr_across_sites, lfp_tfa_avg_tfr_across_sessions
%       diff_condition  - the conditions between which difference has to be
%       calculated, see settings/lfp_tfa_settings_example
%       stat_test       - flag which indicate whether to perform a
%       statistical significance test of differences, set to true while
%       computing difference of averages across site averages of multiple
%       sessions
%		lfp_tfa_cfg     - struct containing the required settings (required
%		only if stat_test = true), see settings/lfp_tfa_settings_example
%           Required fields:
%           fd_rate: Desired false discovery rate
%           fdr_method: method to be used for statistical significance test
%
% OUTPUTS:
%		diff_tfr        - struct containing the LFP time freq spectrogram
%       difference average between different conditions
%
% REQUIRES:
%
% See also lfp_tfa_plot_site_average_tfr,
% lfp_tfa_avg_tfr_across_sites, lfp_tfa_avg_tfr_across_sessions
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

diff_tfr = [];

% defaults
fd_rate = 0.05;
fdr_method = 'pdep';
if nargin < 3
	stat_test = false;
elseif nargin > 3
	stat_test = lfp_tfa_cfg.plot_significant;
	if isfield(lfp_tfa_cfg, 'fd_rate')
		fd_rate = lfp_tfa_cfg.fd_rate;
	end
	if isfield(lfp_tfa_cfg, 'fd_rate')
		fdr_method = lfp_tfa_cfg.fdr_method;
	end
end

if ~isempty([lfp_tfr.cfg_condition])
	conditions = [lfp_tfr.cfg_condition];
else
	return;
end

for i = 1:length(diff_condition)/2
	diff_tfr = struct();
	compare.field = diff_condition{1, 2*(i-1) + 1};
	compare.values = diff_condition{1, 2*(i-1) + 2};
	% check if conditions to compare exist
	if strcmp(compare.field, 'perturbation')
		if sum([compare.values{:}] == unique([conditions.perturbation])) < 2
			continue;
		end
	elseif strcmp(compare.field, 'choice')
		if sum([compare.values{:}] == unique([conditions.choice])) < 2
			continue;
		end
	elseif strcmp(compare.field, 'type_eff')
		if sum(ismember(vertcat(compare.values{:}), ...
				unique([conditions.type; conditions.effector]', 'rows'), 'rows')) < 2
			continue;
		end
	elseif strcmp(compare.field, 'reach_hands')
		if length(unique([conditions.reach_hands])) < 2
			continue
		end
	elseif strcmp(compare.field, 'reach_spaces')
		if length(unique([conditions.reach_spaces])) < 2
			continue
		end
	end
	
	dcn = 0;
	traversed_idx = [];
	for cn = 1:length(lfp_tfr)
		condition_found = false;
		if strcmp(compare.field, 'choice')
			condition_found = lfp_tfr(cn).cfg_condition.choice == compare.values{1};
			
		elseif strcmp(compare.field, 'perturbation')
			condition_found = lfp_tfr(cn).cfg_condition.perturbation == compare.values{1};
			
		elseif strcmp(compare.field, 'type_eff')
			condition_found = lfp_tfr(cn).cfg_condition.type == compare.values{1}(1) ...
				& lfp_tfr(cn).cfg_condition.effector == compare.values{1}(2);
			
		elseif strcmp(compare.field, 'reach_hands')
			condition_found = strcmp(unique(lfp_tfr(cn).cfg_condition.reach_hands), compare.values);
			
		elseif strcmp(compare.field, 'reach_spaces')
			condition_found = strcmp(unique(lfp_tfr(cn).cfg_condition.reach_spaces), compare.values);
		end
		% initially load the pre-injection data structure
		if condition_found
			traversed_idx = [traversed_idx cn];
		else
			continue;
		end
		
		if strcmp(compare.field, 'choice') || strcmp(compare.field, 'perturbation') || strcmp(compare.field, 'type_eff') % the logic is different for hand or space comparisons
			for d = 1:length(lfp_tfr)
				if any(traversed_idx == d), continue; end
				comparison_pair_found = false;
				
				if strcmp(compare.field, 'choice')
					comparison_pair_found = lfp_tfr(d).cfg_condition.type == lfp_tfr(cn).cfg_condition.type ...
						& lfp_tfr(d).cfg_condition.effector == lfp_tfr(cn).cfg_condition.effector ...
						& lfp_tfr(d).cfg_condition.choice == compare.values{2} ...
						& lfp_tfr(d).cfg_condition.perturbation == lfp_tfr(cn).cfg_condition.perturbation;
					
				elseif strcmp(compare.field, 'perturbation')
					comparison_pair_found = lfp_tfr(d).cfg_condition.type == lfp_tfr(cn).cfg_condition.type ...
						& lfp_tfr(d).cfg_condition.effector == lfp_tfr(cn).cfg_condition.effector ...
						& lfp_tfr(d).cfg_condition.choice == lfp_tfr(cn).cfg_condition.choice ...
						& lfp_tfr(d).cfg_condition.perturbation == compare.values{2};
					
				elseif strcmp(compare.field, 'type_eff')
					comparison_pair_found = lfp_tfr(d).cfg_condition.type == compare.values{2}(1) ...
						& lfp_tfr(d).cfg_condition.effector == compare.values{2}(2) ...
						& lfp_tfr(d).cfg_condition.choice == lfp_tfr(cn).cfg_condition.choice ...
						& lfp_tfr(d).cfg_condition.perturbation == lfp_tfr(cn).cfg_condition.perturbation;
				end
				if comparison_pair_found
					
					dcn = dcn + 1;
					%diff_tfr.difference(dcn) = struct();
					% pre injection
					preinj_sync = lfp_tfr(cn);
					if isempty(fieldnames(preinj_sync.hs_tuned_tfs) )
						continue;
					end
					% post injection
					postinj_tfr = lfp_tfr(d);
					if isempty(fieldnames(postinj_tfr.hs_tuned_tfs))
						continue;
					end
					
					if ~isfield(postinj_tfr.hs_tuned_tfs, 'freq') || ...
							~isfield(preinj_sync.hs_tuned_tfs, 'freq')
						continue;
					end
					
					diff_tfr.difference(dcn) = postinj_tfr;
					
					% change the condition label
					diff_tfr.difference(dcn).label = ['( ' postinj_tfr.label, ' - ', ...
						preinj_sync.label ' )'];
					
					diff_tfr.difference(dcn).cfg_condition = postinj_tfr.cfg_condition;
					if strcmp(compare.field, 'choice')
						%                         diff_tfr.difference(dcn).cfg_condition.choice = ['diff' num2str(i)];
						diff_tfr.difference(dcn).cfg_condition.diff = 'choice';
					elseif strcmp(compare.field, 'perturbation')
						%                         diff_tfr.difference(dcn).cfg_condition.perturbation = ['diff' num2str(i)];
						diff_tfr.difference(dcn).cfg_condition.diff = 'perturbation';
					elseif strcmp(compare.field, 'type_eff')
						%                         diff_tfr.difference(dcn).cfg_condition.type_eff = ['diff' num2str(i)];
						diff_tfr.difference(dcn).cfg_condition.diff = 'type_eff';
					end
					
					% loop through handspace tunings
					diff_tfr.difference(dcn).hs_tuned_tfs = postinj_tfr.hs_tuned_tfs;
					for hs = 1:size(postinj_tfr.hs_tuned_tfs, 2)
						for st = 1:size(postinj_tfr.hs_tuned_tfs, 1)
							
							% 20220503sm: make sure to not try to access
							% st or hs values non-existent in preinj_sync.hs_tuned_tfs	
							if ~((st <= size(preinj_sync.hs_tuned_tfs, 1)) && (hs <= size(preinj_sync.hs_tuned_tfs, 2))) 	
								hs_tuned_tfs = postinj_tfr.hs_tuned_tfs;
								for i_hs = 1:size(postinj_tfr.hs_tuned_tfs, 2)
									for i_st = 1:size(postinj_tfr.hs_tuned_tfs, 1)
										
										if (i_st <= size(preinj_sync.hs_tuned_tfs, 1)) && (i_hs <= size(preinj_sync.hs_tuned_tfs, 2))
											hs_tuned_tfs(i_st, i_hs) = preinj_sync.hs_tuned_tfs(i_st, i_hs);
										end	
									end
									preinj_sync.hs_tuned_tfs = hs_tuned_tfs;
								end
							end
							
							
							if (st <= size(preinj_sync.hs_tuned_tfs, 1)) && (hs <= size(preinj_sync.hs_tuned_tfs, 2)) && ...
									isfield(preinj_sync.hs_tuned_tfs(st, hs).freq, 'powspctrm') && ...
									isfield(postinj_tfr.hs_tuned_tfs(st, hs).freq, 'powspctrm') && ...
									~isempty(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm) && ...
									~isempty(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm)
								ntimebins = min([size(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm, 3), ...
									size(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm, 3)]);
								diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.time = ...
									postinj_tfr.hs_tuned_tfs(st, hs).freq.time(1:ntimebins);
								% calculate the difference
								if isfield(postinj_tfr.hs_tuned_tfs(st, hs), 'ntrials')
									diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
										nanmean(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1) - ...
										nanmean(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1);
									diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).ntrials = ...
										[];
								else
									diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
										postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins) - ...
										preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins);
									% statistical significance test
									if stat_test == true
										% paired ttest
										[~, p] = ttest(...
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm);
										% multiple comparison correction
										if strcmp(lfp_tfa_cfg.correction_method,'FDR');
											[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
												fdr_method, 'yes');
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.h = h;
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.crit_p = crit_p;
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.adj_ci_cvrg = adj_ci_cvrg;
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.adj_p = adj_p;
										elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
											corrected_p_val_sig = 0.05/(size(preinj_sync.hs_tuned_tfs(st, hs).freq.time,2)*size(preinj_sync.hs_tuned_tfs(st, hs).freq.freq,2));
											diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.h = p < corrected_p_val_sig;
										end
									end
								end
								
							else
								diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = [];
								diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.time = [];
								diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.freq = [];
							end
						end
					end
				else
					continue;
				end
			end
			
		elseif strcmp(compare.field, 'reach_hands') || strcmp(compare.field, 'reach_spaces') % need to average same hand or same space for each lfp_tfr(d), then do the difference between the 2 averages
			
			
			for d = 1:length(lfp_tfr)
				diff_tfr.difference(d) = lfp_tfr(d);
				if strcmp(compare.field, 'reach_hands')
					diff_tfr.difference(d).hs_tuned_tfs = lfp_tfr(d).hs_tuned_tfs;
					diff_tfr.difference(d).label = [lfp_tfr(d).label, 'CH - IH'];
					%                     diff_tfr.difference(d).cfg_condition.reach_hands = ['diff' num2str(i)];
					diff_tfr.difference(d).cfg_condition = lfp_tfr(d).cfg_condition;
					diff_tfr.difference(d).cfg_condition.diff = 'hands';
					
					for st = 1:size(lfp_tfr(d).hs_tuned_tfs, 1) %st is windows here
						%calculate difference between same space,
						%opposite hands condition
						ntimebins_CS = min([size(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm, 3)...
							size(lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm, 3)]);
						ntimebins_IS = min([size(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm, 3)...
							size(lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm, 3)]);
						
						if isfield(lfp_tfr(d).hs_tuned_tfs, 'ntrials')
							diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = nanmean(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_CS))...
								- nanmean(lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm(:,:,1:ntimebins_CS));
							diff_tfr.difference(d).hs_tuned_tfs(st, 2).freq.powspctrm = nanmean(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_IS))...
								- nanmean(lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm(:,:,1:ntimebins_IS));
							
							
						else
							
							
							
							diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_CS)...
								- lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm(:,:,1:ntimebins_CS);
							diff_tfr.difference(d).hs_tuned_tfs(st, 2).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_IS)...
								- lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm(:,:,1:ntimebins_IS);
							
						end
						%change hand space label
						diff_tfr.difference(d).hs_tuned_tfs(st, 1).hs_label = {'CHCS - IHCS'};
						diff_tfr.difference(d).hs_tuned_tfs(st, 2).hs_label = {'CHIS - IHIS'};
						%empty other hand space condition since not needed
						%here
						diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.powspctrm = [];
						diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.time = [];
						diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.freq = [];
						diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.powspctrm = [];
						diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.time = [];
						diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.freq = [];
						
						
						if stat_test == true
							for hs = 1:2
								% paired ttest
								[~, p] = ttest(...
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.powspctrm);
								% multiple comparison correction
								if strcmp(lfp_tfa_cfg.correction_method,'FDR');
									[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
										fdr_method, 'yes');
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.h = h;
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.crit_p = crit_p;
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.adj_ci_cvrg = adj_ci_cvrg;
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.adj_p = adj_p;
								elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
									corrected_p_val_sig = 0.05/size(lfp_tfr(d).hs_tuned_tfs(st,hs).freq.powspctrm,1);
									diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.h = p < corrected_p_val_sig;
								end
							end
						end
					end
					
				elseif strcmp(compare.field, 'reach_spaces')
					
					diff_tfr.difference(d).hs_tuned_tfs = lfp_tfr(d).hs_tuned_tfs;
					diff_tfr.difference(d).label = [lfp_tfr(d).label, 'CS - IS'];
					%                     diff_tfr.difference(d).cfg_condition.reach_spaces = ['diff' num2str(i)];
					diff_tfr.difference(d).cfg_condition = lfp_tfr(d).cfg_condition;
					diff_tfr.difference(d).cfg_condition.diff = 'spaces';
					
					
					for st = 1:size(lfp_tfr(d).hs_tuned_tfs, 1) %st is windows here
						%calculate difference between same hands,
						%opposite space condition
						if ismember(lfp_tfr(d).cfg_condition.effector,[1 2 3 4 6])
							
							
							if size(lfp_tfr(d).hs_tuned_tfs, 2) == 1
								continue
							end
							
							% try to figure out which reach_hands were in
							% use
							unique_reach_hand_list = unique(lfp_tfr(d).cfg_condition.reach_hands);
							n_reach_hands = length(unique_reach_hand_list);
							
							
							% check for empty handpaces and fill with
							% NaNed filled handspace data
							empty_handspace_list = zeros([1, size(lfp_tfr(d).hs_tuned_tfs, 2)]);
							for i_hand_reach_combos = 1 : size(lfp_tfr(d).hs_tuned_tfs, 2)
								if isempty(lfp_tfr(d).hs_tuned_tfs(st, i_hand_reach_combos).freq)
									empty_handspace_list(i_hand_reach_combos) = 1;
								end
							end
							empty_handspace_idx = find(empty_handspace_list);
							
							hs_label_list = {'CH CS', 'CH CI', 'IH CS', 'IH IS'};
							
							nonempty_handspace_idx = find(empty_handspace_list == 0);
							first_non_empty_handspace_idx = nonempty_handspace_idx(1);
							if ~isempty(empty_handspace_idx)
								for i_empty_handspace = 1 : length(empty_handspace_idx)
									lfp_tfr(d).hs_tuned_tfs(st, empty_handspace_idx(i_empty_handspace)) = lfp_tfr(d).hs_tuned_tfs(st, first_non_empty_handspace_idx);
									lfp_tfr(d).hs_tuned_tfs(st, empty_handspace_idx(i_empty_handspace)).hs_label = hs_label_list(empty_handspace_idx(i_empty_handspace));
									lfp_tfr(d).hs_tuned_tfs(st, empty_handspace_idx(i_empty_handspace)).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st, empty_handspace_idx(i_empty_handspace)).freq.powspctrm * NaN;
								end
							end
							
							if (n_reach_hands == 1)
								disp(['Found only read_hand ', unique_reach_hand_list{1}, ' duplicating to ''fake'' the missing reach_hand for now.']);
								
								% 								if size(lfp_tfr(d).hs_tuned_tfs, 2) > 2
								% 									% this should not happen, if only
% 									% either CH or IH trials exist we
% 									% assume their CS IS will be in 1 and 2
% 									% respectively
% 									keyboard
% 								end								

								% make sure to create the artifical data so
								% the difference operation will work
								lfp_tfr(d).hs_tuned_tfs(st, 3) =  lfp_tfr(d).hs_tuned_tfs(st, 1);
								lfp_tfr(d).hs_tuned_tfs(st, 4) = lfp_tfr(d).hs_tuned_tfs(st, 2);
								% adjust labels and NaN out the artificial
								% data so the plots are easy to recognize
								% as invalid
								

								
								switch unique_reach_hand_list{1}
									case 'C'
										lfp_tfr(d).hs_tuned_tfs(st, 3).hs_label = {'IH CS'};
										lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm * NaN;
								
										lfp_tfr(d).hs_tuned_tfs(st, 4).hs_label = {'IH IS'};
										lfp_tfr(d).hs_tuned_tfs(st, 4).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm * NaN;
									case 'I'
										lfp_tfr(d).hs_tuned_tfs(st, 1).hs_label = {'CH CS'};
										lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm * NaN;
										
										lfp_tfr(d).hs_tuned_tfs(st, 1).hs_label = {'CH IS'};
										lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm * NaN;
										
								end
							end
							
							ntimebins_CH = min([size(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm, 3)...
								size(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm, 3)]);
							ntimebins_IH = min([size(lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm, 3)...
								size(lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm, 3)]);
							
							
							if isfield(lfp_tfr(d).hs_tuned_tfs, 'ntrials')
								diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = nanmean(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_CH))...
									- nanmean(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_CH));
								diff_tfr.difference(d).hs_tuned_tfs(st, 2).freq.powspctrm = nanmean(lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm(:,:,1: ntimebins_IH))...
									- nanmean(lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm(:,:,1: ntimebins_IH));
								
							else
								diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_CH)...
									- lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_CH);
								diff_tfr.difference(d).hs_tuned_tfs(st, 2).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,3).freq.powspctrm(:,:,1: ntimebins_IH)...
									- lfp_tfr(d).hs_tuned_tfs(st,4).freq.powspctrm(:,:,1: ntimebins_IH);
							end
							%change hand space label
							diff_tfr.difference(d).hs_tuned_tfs(st, 1).hs_label = {'CHCS - CHIS'};
							diff_tfr.difference(d).hs_tuned_tfs(st, 2).hs_label = {'IHCS - IHIS'};
							%empty other hand space condition since not needed
							%here
							diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.powspctrm = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.time = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,3).freq.freq = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.powspctrm = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.time = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,4).freq.freq = [];
							
							
							if stat_test == true
								for hs = 1:2
									if (sum(isnan(diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.powspctrm(:)))) > 0 || isempty(diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.powspctrm)
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.h = 0;
										continue
									end
									% paired ttest
									[~, p] = ttest(...
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.powspctrm);
									% multiple comparison correction
									if strcmp(lfp_tfa_cfg.correction_method,'FDR');
										[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
											fdr_method, 'yes');
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.h = h;
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.crit_p = crit_p;
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.adj_ci_cvrg = adj_ci_cvrg;
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.adj_p = adj_p;
									elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
										corrected_p_val_sig = 0.05/size(lfp_tfr(d).hs_tuned_tfs(st,hs).freq.powspctrm,1);
										diff_tfr.difference(d).hs_tuned_tfs(st, hs).freq.stat_test.h = p < corrected_p_val_sig;
									end
									
								end
							end
						elseif lfp_tfr(d).cfg_condition.effector == 0
							ntimebins_space = min([size(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm, 3)...
								size(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm, 3)]);
							
							if isfield(lfp_tfr(d).hs_tuned_tfs, 'ntrials')
								diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = nanmean(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_space),1)...
									- nanmean(lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_space),1);
								
								
							else
								diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm = lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm(:,:,1:ntimebins_space)...
									- lfp_tfr(d).hs_tuned_tfs(st,2).freq.powspctrm(:,:,1:ntimebins_space);
								
							end
							%change hand space label
							diff_tfr.difference(d).hs_tuned_tfs(st, 1).hs_label = {'anyH_CS - anyH_IS'};
							%empty other hand space condition since not needed
							%here
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.powspctrm = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.time = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.freq = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.powspctrm = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.time = [];
							diff_tfr.difference(d).hs_tuned_tfs(st,2).freq.freq = [];
							
							if stat_test == true
								
								% paired ttest
								[~, p] = ttest(...
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.powspctrm);
								% multiple comparison correction
								if strcmp(lfp_tfa_cfg.correction_method,'FDR');
									[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
										fdr_method, 'yes');
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.stat_test.h = h;
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.stat_test.crit_p = crit_p;
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.stat_test.adj_ci_cvrg = adj_ci_cvrg;
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.stat_test.adj_p = adj_p;
								elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
									corrected_p_val_sig = 0.05/size(lfp_tfr(d).hs_tuned_tfs(st,1).freq.powspctrm,1);
									diff_tfr.difference(d).hs_tuned_tfs(st, 1).freq.stat_test.h = p < corrected_p_val_sig;
								end
								
								
							end
							
							
						end
					end
					
				end
				
				
			end
			
			
			
		end
	end
	if isfield(diff_tfr, 'difference')
		lfp_tfr = diff_tfr.difference;
	end
end
% generate return variable
if isfield(diff_tfr, 'difference')
	diff_tfr = diff_tfr.difference;
else
	diff_tfr = [];
end
end

