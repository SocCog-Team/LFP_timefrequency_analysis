function [ cmp_conditions ] = lfp_tfa_compare_conditions( lfp_tfa_cfg, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    task_types = lfp_tfa_cfg.compare.types;
    effectors = lfp_tfa_cfg.compare.effectors;
    targets = lfp_tfa_cfg.compare.targets;
    if nargin > 1
        targets = varargin{1};
    end
    %unique([states_lfp.recorded_hemispace]);
    if lfp_tfa_cfg.compare.choice_trials
        choices = unique([states_lfp(1).trials.choice_trial]);
    else
        choices = lfp_tfa_cfg.compare.choice_trials;
    end
    perturbations = lfp_tfa_cfg.compare.perturbations;
    perturbation_groups = lfp_tfa_cfg.compare.perturbation_groups;
    
    hands = lfp_tfa_cfg.compare.reach_hands;
    spaces = lfp_tfa_cfg.compare.reach_spaces;  
    % assign hand space labels
    hs_labels = cell(1, length(hands)*length(spaces));
    reach_hands = cell(1, length(hands)*length(spaces));
    reach_spaces = cell(1, length(hands)*length(spaces));
    for h = 1:length(hands)
        if strcmp(hands{h},'R') || strcmp(hands{h},'L')
            if strcmp(hands{h},lfp_tfa_cfg.ref_hemisphere)
                hand_label = 'IH';
            else
                hand_label = 'CH';
            end
        else
            hand_label = [hands{h}, 'H'];
        end
        for s = 1:length(spaces)
            if strcmp(spaces{s},'R') || strcmp(spaces{s},'L')
                if strcmp(spaces{s},lfp_tfa_cfg.ref_hemisphere)
                    space_label = 'IS';
                else
                    space_label = 'CS';
                end
            else
                space_label = [spaces{s}, 'S'];
            end
            hs_idx = (h-1)*length(spaces) + s;
            reach_hands{hs_idx} = hands{h};
            reach_spaces{hs_idx} = spaces{s};
            hs_labels{hs_idx} = [hand_label ' ' space_label];
        end
    end
    
    % create conditions
    cmp_conditions = struct();
    
    i = 0;
    %for target = targets
        %target_label = target{1};
        for type = task_types
            type_label = ['Type_' num2str(type)];
            for eff = effectors
                eff_label = ['Eff_' num2str(eff)];            
                for ch = choices
                    if ch == 0
                        ch_label = 'Instr';
                    elseif choice == 1
                        ch_label = 'Choice';
                    else
                        ch_label = [];
                    end
                    for p = 1:length(perturbations)
                        if perturbations(p) == 0
                            p_label = 'Pre';
                        elseif perturbations(p) == 1
                            p_label = 'Post';
                        else
                            p_label = [];
                        end 
                        
                        i = i + 1;
                        condition_label = [type_label, '_', eff_label, '_', ...
                            ch_label, '_', p_label];
                        cmp_conditions(i).type = type;
                        cmp_conditions(i).effector = eff;
                        %cmp_conditions(i).target = target{1};
                        cmp_conditions(i).choice = ch;
                        cmp_conditions(i).perturbation = perturbations(p);
                        cmp_conditions(i).perturbation_group = perturbation_groups(p);
                        cmp_conditions(i).hs_labels = hs_labels;
                        cmp_conditions(i).reach_hands = reach_hands;
                        cmp_conditions(i).reach_spaces = reach_spaces;
                        cmp_conditions(i).label = condition_label;
%                             end
%                         end
                    end
                end
            end
        end
    %end
                        
                        %     for rec_hem = recorded_hemispace
%         for c = choice
%             for b = blocks
%                 i = i + 1;
%                 cfg_conditions(i).recorded_hemispace = rec_hem;
%                 cfg_conditions(i).choice = c;
%                 cfg_conditions(i).block = b;
%                 cfg_conditions(i).perturbation = ...
%                     unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
%                 cond_label = [];
%                 if cfg_conditions(i).recorded_hemispace == 'L'
%                     cond_label = [cond_label 'Left_hemisphere_'];
%                 else
%                     cond_label = [cond_label 'Right_hemisphere_'];
%                 end
%                 if cfg_conditions(i).choice == 0
%                     cond_label = [cond_label 'Instructed_'];
%                 else
%                     cond_label = [cond_label 'Choice_'];
%                 end
%                 if cfg_conditions(i).perturbation == 0
%                     cond_label = [cond_label 'Control_'];
%                 else
%                     cond_label = [cond_label 'Inactivation_'];
%                 end
%                 cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
%                 cfg_conditions(i).label = cond_label;
%                                                 
%             end
%         end
%     end
    
%     i = 0;
% 
%         
%     
%     
%     if isfield(lfp_tfa_cfg, 'add_conditions')
%         for c = 1:length(lfp_tfa_cfg.add_conditions)
%             if ~isempty(lfp_tfa_cfg.add_conditions(c))
%                 if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
%                     for rec_hem = recorded_hemispace        
%                         for ch = choice
%                             i = i + 1;
%                             if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'control')
%                                 cmp_conditions(i).block = blocks(perturbation == 0);
%                                 cmp_conditions(i).perturbation = 0;
%                             elseif strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
%                                 cmp_conditions(i).block = blocks(perturbation ~= 0);
%                                 cmp_conditions(i).perturbation = 1;
%                             else
%                                 cmp_conditions(i).block = blocks(lfp_tfa_cfg.add_conditions(c).blocks);
%                                 cmp_conditions(i).perturbation = sign(perturbation(blocks == blocks(lfp_tfa_cfg.add_conditions(c).blocks(1))));
%                                 
%                             end                    
%                             cmp_conditions(i).choice = ch;
%                             if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
%                                 cmp_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
%                             end
%                             cmp_conditions(i).recorded_hemispace = rec_hem;
%                             cond_label = [];
%                             if cmp_conditions(i).recorded_hemispace == 'L'
%                                 cond_label = [cond_label 'Left hemisphere_'];
%                             else
%                                 cond_label = [cond_label 'Right hemisphere_'];
%                             end
%                             if cmp_conditions(i).choice == 0
%                                 cond_label = [cond_label 'Instructed Trials_'];
%                             else
%                                 cond_label = [cond_label 'Choice Trials_'];
%                             end
%                             if cmp_conditions(i).perturbation == 0
%                                 cond_label = [cond_label 'Pre-Injection_'];
%                             else
%                                 cond_label = [cond_label 'Post-Injection_'];
%                             end
%                             cond_label = [cond_label, 'Block_', num2str(cmp_conditions(i).block)];
%                             cmp_conditions(i).label = cond_label;
%                             
%                         end
%                     end
%                 end
%             end
%         end
%     end

end

