function sites_avg = lfp_tfa_avg_pow_across_sites(lfp_pow, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_pow_across_sites  - Condition-based LFP power spectrum
%(Power vs. Frequency) grand average across many site averages from
%multiple sessions
%
% USAGE:
%	sites_avg = lfp_tfa_avg_pow_across_sites(lfp_evoked, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based average LFP
%       power spectrum for	multiple sites, output of lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               1. session.sites - 1xM struct whose each element in turn is 
%               a 1xN structs containing condition-based average LFP power 
%               spectrum for N sites recorded in one session (M = number of
%               sessions analysed)
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               settings/lfp_tfa_settings_example and 
%               lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are 
%               saved. Results will be saved under 
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sites/LFP_Power']
%               3. compare.targets     - target areas to compare, see 
%               settings/lfp_tfa_settings_example
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example
% OUTPUTS:
%		sites_avg       - structure containing condition-based LFP Power
%		spectrum response averaged across multiple sites
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_psd
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings, 
% lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_psd
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

    % results folder
    if nargin < 3
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_Power');
    else
        results_fldr = varargin{1};
    end
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average Evoked LFP across sites
    sites_avg = struct();
    
    for t = 1:length(lfp_tfa_cfg.compare.targets)
        sites_avg(t).target = lfp_tfa_cfg.compare.targets{t};
        for cn = 1:length(lfp_tfa_cfg.conditions)
            fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
            sites_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            %sites_avg(t).condition(cn).avg_across_sessions = struct();
            % initialize number of sites for each handspace
            % label
            for ep = 1:size(lfp_pow.session(1).sites(1).condition(cn).hs_tuned_power, 1)
                for hs = 1:size(lfp_pow.session(1).sites(1).condition(cn).hs_tuned_power, 2)
                    sites_avg(t).condition(cn).hs_tuned_power(ep, hs).nsites = 0;
                    sites_avg(t).condition(cn).hs_tuned_power(ep, hs).psd = [];
                end
            end
            nsites = 0;
            for i = 1:length(lfp_pow.session)  
                for j = 1:length(lfp_pow.session(i).sites)
                    if ~strcmp(lfp_pow.session(i).sites(j).target, lfp_tfa_cfg.compare.targets{t})
                        continue;
                    end
                    if ~lfp_pow.session(i).sites(j).use_for_avg
                        continue;
                    end
                    if ~isempty(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power) && ... 
                        isfield(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 'mean')
                        nsites = nsites + 1;   
                        for ep = 1:size(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 1)
                            for hs = 1:size(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 2)
                                if isfield(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs), 'mean') ...
                                        && ~isempty(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean)
                                    sites_avg(t).condition(cn).hs_tuned_power(ep, hs).nsites = ...
                                        sites_avg(t).condition(cn).hs_tuned_power(ep, hs).nsites + 1;
                                    if sites_avg(t).condition(cn).hs_tuned_power(ep, hs).nsites == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).freq ...
                                        = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq;
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).hs_label ...
                                            = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).hs_label;
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).epoch_name ...
                                            = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).epoch_name;
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).psd ...
                                            = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean;
                                    else
                                        nfreqbins = length(sites_avg(t).condition(cn).hs_tuned_power(ep, hs).freq);
                                        % average same number of time bins
                                        if nfreqbins > length(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq)
                                            nfreqbins = length(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq);
                                        end
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).freq = ...
                                            lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq(1:nfreqbins);
                                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).psd ...
                                            = [lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqbins); ...
                                            sites_avg(t).condition(cn).hs_tuned_power(ep,hs).psd(:,1:nfreqbins)];                              
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % compute average
            for ep = 1:size(sites_avg(t).condition(cn).hs_tuned_power, 1)
                for hs = 1:size(sites_avg(t).condition(cn).hs_tuned_power, 2)
                    if isfield(sites_avg(t).condition(cn).hs_tuned_power(ep,hs), 'psd')
                        sites_avg(t).condition(cn).hs_tuned_power(ep,hs).dimord = 'nsites_freq';                                
                        [sites_avg(t).condition(cn).hs_tuned_power(ep,hs).mean, ...
                            sites_avg(t).condition(cn).hs_tuned_power(ep,hs).error] = ...
                            (lfp_tfa_compute_statistics(...
                            sites_avg(t).condition(cn).hs_tuned_power(ep,hs).psd, lfp_tfa_cfg.error_measure));
                    end
                end
            end


            if ~isempty(sites_avg(t).condition(cn).hs_tuned_power)
                if isfield(sites_avg(t).condition(cn).hs_tuned_power,... 
                        'mean')
                    plottitle = [lfp_tfa_cfg.compare.targets{t}, ...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        lfp_tfa_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    sprintf('LFP_Power_%s_%s', lfp_tfa_cfg.compare.targets{t}, lfp_tfa_cfg.conditions(cn).label ));
                    lfp_tfa_plot_hs_tuned_psd_2(sites_avg(t).condition(cn).hs_tuned_power, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end
            
        end
    end
    % save session average power spectrum
    save(fullfile(results_fldr, 'LFP_Power_sites_avg.mat'), 'sites_avg');
    close all;
end