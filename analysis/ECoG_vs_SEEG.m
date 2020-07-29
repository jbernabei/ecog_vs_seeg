clear all

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field,... 
    id_field,implant_field, outcome_field, resect_field, roi_field,...
    target_field, therapy_field] = set_up_workspace(iEEG_atlas_path);

%% 1. Patient-specific base network renderings

good_ecog_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good') & strcmp(implant_field,'ECoG'));
good_seeg_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good') & strcmp(implant_field,'SEEG'));

ecog_patients = all_patients(good_ecog_patient_indices);
seeg_patients = all_patients(good_seeg_patient_indices);

rendering_filename = 'base_rendering.png';

%create_rendering(rendering_filename, {ecog_patients.coords}, {ecog_patients.conn}, {ecog_patients.patientID});
%create_rendering(rendering_filename, {seeg_patients.coords}, {seeg_patients.conn}, {seeg_patients.patientID});

%% 2A. Anatomical targets, broken down by implant type and target (Frontal, Temporal, other)

ecog_temporal_indices = find([hasData_field{:}] & strcmp(outcome_field,'good')...
    & strcmp(implant_field,'ECoG') & (strcmp(target_field,'Temporal') | strcmp(target_field,'MTL')));
ecog_frontal_indices = find([hasData_field{:}] & strcmp(outcome_field,'good')...
    & strcmp(implant_field,'ECoG') & (strcmp(target_field,'Frontal') | strcmp(target_field,'MFL') | strcmp(target_field,'FP')));

seeg_temporal_indices = find([hasData_field{:}] & strcmp(outcome_field,'good')...
    & strcmp(implant_field,'SEEG') & (strcmp(target_field,'Temporal') | strcmp(target_field,'MTL')));
seeg_frontal_indices = find([hasData_field{:}] & strcmp(outcome_field,'good')...
    & strcmp(implant_field,'SEEG') & (strcmp(target_field,'Frontal') | strcmp(target_field,'MFL') | strcmp(target_field,'FP')));

[ECoG_top_regions, ECoG_plot_data] = rank_anatomical_targets({ecog_patients.patientID},{ecog_patients.roi},...
    {ecog_patients.laterality}, all_inds, all_locs,'ECoG','Temporal','plot1.fig');

[SEEG_top_regions, SEEG_plot_data] = rank_anatomical_targets({seeg_patients.patientID},{seeg_patients.roi},...
    {seeg_patients.laterality}, all_inds, all_locs,'ECoG','Frontal','plot2.fig');

% this needs some small fixes
figure(1);clf;
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
subplot(1,2,1)
bar(ET_plot_data(1:15))
title('Anatomical distribution of ECoG patients');
set(gca,'xtick',(1:15),'xticklabel',ECoG_top_regions(1:15));
xtickangle(45);
subplot(1,2,2)
bar(EF_plot_data(1:15))
title('Anatomical distribution of SEEG patients');
set(gca,'xtick',(1:15),'xticklabel',SEEG_top_regions(1:15));
xtickangle(45);

%% 2B. Interelectrode distance comparison

[ecog_distances_all, ecog_distances_pt] = compute_interelectrode_distances({ecog_patients.coords});
[seeg_distances_all, seeg_distances_pt] = compute_interelectrode_distances({seeg_patients.coords});

figure(2);clf
subplot(1,2,1)
histogram(ecog_distances_all)
subplot(1,2,2)
histogram(seeg_distances_all)

%% 2C. Connectivity analysis

% histogram of raw values / edge variability for each node -> lets do HG
% break down each histogram / bar by region (intra ROI / inter ROI)

[ecog_conn_all, ecog_conn_pt] = analyze_connectivity({ecog_patients.conn});
[seeg_conn_all, seeg_conn_pt] = analyze_connectivity({seeg_patients.conn});

figure(2);clf
subplot(5,1,1)
hold on
histogram(ecog_conn_all(:,1))
histogram(seeg_conn_all(:,1))
subplot(5,1,2)
hold on
histogram(ecog_conn_all(:,2))
histogram(seeg_conn_all(:,2))
subplot(5,1,3)
hold on
histogram(ecog_conn_all(:,3))
histogram(seeg_conn_all(:,3))
subplot(5,1,4)
hold on
histogram(ecog_conn_all(:,4))
histogram(seeg_conn_all(:,4))
subplot(5,1,5)
hold on
histogram(ecog_conn_all(:,5))
histogram(seeg_conn_all(:,5))


%% 3A. Small worldness (or other global metric)

a = 0;
for pt = good_ecog_patient_indices
    a = a+1;
    for f = 1:5
        pt_adj = all_patients(pt).conn(f).data;
        [SWP_ecog,delta_C_ecog,delta_L_ecog] = small_world_propensity(pt_adj);
        sw_value_ecog(a,f) = SWP_ecog;
    end
end

a = 0;
for pt = good_seeg_patient_indices
    a = a+1;
    for f = 1:5
        pt_adj = all_patients(pt).conn(f).data; 
        [SWP_seeg,delta_C_seeg,delta_L_seeg] = small_world_propensity(pt_adj);
        sw_value_seeg(a,f) = SWP_seeg;
    end
end

figure(3);clf;
subplot(1,2,1)
boxplot(sw_value_ecog)
ylabel('Small-world propensity')
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylim([0.25 1])
title('Small-worldness of ECoG networks')
subplot(1,2,2)
boxplot(sw_value_seeg)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Small-world propensity')
title('Small-worldness of SEEG networks')
ylim([0.25 1])

%% 3B. Modularity analysis -> compute modularity & participation coefficient

[ecog_pt_modules, ecog_pt_q_vals, ecog_pt_pc] = compute_modularity({ecog_patients.conn});
[seeg_pt_modules, seeg_pt_q_vals, seeg_pt_pc] = compute_modularity({seeg_patients.conn});

%% 3C. Comparison of modules to Yeo systems & own electrode groups

% have to get 'community' labels for electrodes and Yeo systems

% must do a network similarity function

%% 4A. Node strength of the SOZ/EZ

[ecog_str_non_res, ecog_str_res, ecog_str] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'node_strength');
[seeg_str_non_res, seeg_str_res, seeg_str] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'node_strength');

figure(4);clf;
subplot(2,2,1)
boxplot(ecog_str_non_res)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Z score node strength')
title('Non-resected node strength ECoG')
ylim([-1 3])
subplot(2,2,2)
boxplot(ecog_str_res)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Z score node strength')
title('Resected node strength ECoG')
ylim([-1 3])
subplot(2,2,3)
boxplot(seeg_str_non_res)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Z score node strength')
title('Non resected node strength SEEG')
ylim([-1 3])
subplot(2,2,4)
boxplot(seeg_str_res)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Z score node strength')
title('Resected node strength SEEG')
ylim([-1 3])

%% 4B. Betweenness centrality of the SOZ/EZ

[ecog_btw_non_res, ecog_btw_res, ecog_btw] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'betweenness_centrality');
[seeg_btw_non_res, seeg_btw_res, seeg_btw] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'betweenness_centrality');

figure(5);clf;
subplot(2,2,1)
boxplot(ecog_btw_non_res)
set(gca,'xtick',(1:5),'xticklabel',{'Broadband CC','Alpha/Theta','Beta','Low-Gamma','High-Gamma'});
ylabel('Z score betweenness centrality')
title('Non-resected betweenness centrality ECoG')
ylim([-1 3])
subplot(2,2,2)
boxplot(ecog_btw_res)
ylabel('Z score betweenness centrality')
title('Resected betweenness centrality ECoG')
ylim([-1 3])
subplot(2,2,3)
boxplot(seeg_btw_non_res)
ylabel('Z score betweenness centrality')
title('Non-resected betweenness centrality SEEG')
ylim([-1 3])
subplot(2,2,4)
boxplot(seeg_btw_res)
ylabel('Z score betweenness centrality')
title('Resected betweenness centrality SEEG')
ylim([-1 3])
%% 4C. Control centrality of the SOZ/EZ

[ecog_cc_non_res, ecog_cc_res, ecog_cc] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'control_centrality');
[seeg_cc_non_res, seeg_cc_res, seeg_cc] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'control_centrality');

figure(5);clf;
subplot(2,2,1)
boxplot(ecog_cc_non_res)
ylabel('Z score control centrality')
title('Non-resected control centrality ECoG')
ylim([-1 3])
subplot(2,2,2)
boxplot(ecog_cc_res)
ylabel('Z score control centrality')
title('Resected control centrality ECoG')
ylim([-1 3])
subplot(2,2,3)
boxplot(seeg_cc_non_res)
ylabel('Z score control centrality')
title('Non-resected control centrality SEEG')
ylim([-1 3])
subplot(2,2,4)
boxplot(seeg_cc_res)
ylabel('Z score control centrality')
title('Resected control centrality SEEG')
ylim([-1 3])

%% 5A. Localization based on top metric -> All nodes

% need 2 localization methods -> use the localize_EZ_II method -> adapt
[targets_ecog_raw, dice_ecog_raw] = localize_EZ({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.resect}, 'node_strength');
[targets_seeg_raw, dice_seeg_raw] = localize_EZ({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.resect}, 'node_strength');

%% 5B. Localization based on top metric -> No WM

% generate new data
[new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, 'no_WM');
[new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, 'no_WM');

[targets_ecog_raw, dice_ecog_raw] = localize_EZ(adj_matrices, mni_coordinates, resected_elecs, metric_type);
[targets_seeg_raw, dice_seeg_raw] = localize_EZ(adj_matrices, mni_coordinates, resected_elecs, metric_type);

%% 5C. Localization based on top metric -> ROI reduced

% generate new data
[new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, 'min_ROI');
[new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, 'min_ROI');

[targets_ecog_raw, dice_ecog_raw] = localize_EZ(adj_matrices, mni_coordinates, resected_elecs, metric_type);
[targets_seeg_raw, dice_seeg_raw] = localize_EZ(adj_matrices, mni_coordinates, resected_elecs, metric_type);