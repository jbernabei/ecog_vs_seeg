clear all

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field,... 
    id_field,implant_field, outcome_field, resect_field, roi_field,...
    target_field, therapy_field] = set_up_workspace(iEEG_atlas_path);

% set up colors
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];

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

[ECoG_top_regions, ECoG_plot_data, ECoG_total_elecs] = rank_anatomical_targets({ecog_patients.patientID},{ecog_patients.roi},...
    {ecog_patients.laterality}, all_inds, all_locs,'ECoG','Temporal','plot1.fig');

[SEEG_top_regions, SEEG_plot_data, SEEG_total_elecs] = rank_anatomical_targets({seeg_patients.patientID},{seeg_patients.roi},...
    {seeg_patients.laterality}, all_inds, all_locs,'ECoG','Frontal','plot2.fig');

% this needs some small fixes
figure(1);clf;
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
subplot(1,2,1)
plotdata1 = ECoG_plot_data(1:15)./ECoG_total_elecs;
bh1 = bar(1:numel(plotdata1),diag(plotdata1),'stacked','FaceColor', color1);
bh1(12).FaceColor = color2; % contralateral coloration
bh1(13).FaceColor = color2; % contralateral coloration
ylabel('Fraction of total electrodes')
title('Anatomical distribution of ECoG patients');
set(gca,'xtick',(1:15),'xticklabel',ECoG_top_regions(1:15));
ylim([0, 0.11])
xtickangle(45);
subplot(1,2,2)
plotdata2 = SEEG_plot_data(1:15)./SEEG_total_elecs;
bh2 = bar(1:numel(plotdata2),diag(plotdata2),'stacked','FaceColor', color1);
bh2(3).FaceColor = color2; % contralateral coloration
bh2(6).FaceColor = color2; % contralateral coloration
bh2(12).FaceColor = color2; % contralateral coloration
bh2(14).FaceColor = color2; % contralateral coloration
title('Anatomical distribution of SEEG patients');
set(gca,'xtick',(1:15),'xticklabel',SEEG_top_regions(1:15));
ylabel('Fraction of total electrodes')
xtickangle(45);
ylim([0, 0.11])

%% 2B. Interelectrode distance comparison

[ecog_distances_all, ecog_distances_pt] = compute_interelectrode_distances({ecog_patients.coords});
[seeg_distances_all, seeg_distances_pt] = compute_interelectrode_distances({seeg_patients.coords});

figure(2);clf
subplot(1,2,1)
histogram(ecog_distances_all,'Normalization','probability','FaceColor',color1)
xlim([0 200])
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG internodal distances')
subplot(1,2,2)
histogram(seeg_distances_all,'Normalization','probability','FaceColor',color1)
xlim([0 200])
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('SEEG internodal distances')

% calculate means
m1 = mean(ecog_distances_all)
m2 = mean(seeg_distances_all)

% calculate skewness
s1 = skewness(ecog_distances_all)
s2 = skewness(seeg_distances_all)

figure(3);clf;
hold on
histogram(ecog_distances_all,'Normalization','probability','FaceColor',color1)
histogram(seeg_distances_all,'Normalization','probability','FaceColor',color2)
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG versus SEEG internodal distances')
legend('ECoG','SEEG','Location','NorthEast')
hold off
%% 2C. Connectivity analysis

% histogram of raw values / edge variability for each node -> lets do HG
% break down each histogram / bar by region (intra ROI / inter ROI)

[ecog_conn_all, ecog_conn_pt] = analyze_connectivity({ecog_patients.conn});
[seeg_conn_all, seeg_conn_pt] = analyze_connectivity({seeg_patients.conn});

figure(2);clf
subplot(2,1,1)
hold on
histogram(seeg_conn_all(:,1),'Normalization','probability','FaceColor', color2)
histogram(ecog_conn_all(:,1),'Normalization','probability','FaceColor', color1)
legend('SEEG','ECoG','Location','NorthEast')
xlim([0.1,0.6])
title('Broadband cross-correlation edge weight distribution')
ylabel('Frequency')
subplot(2,1,2)
hold on
histogram(seeg_conn_all(:,3),'Normalization','probability','FaceColor', color2)
histogram(ecog_conn_all(:,3),'Normalization','probability','FaceColor', color1)
xlim([0.1,0.6])
title('Beta coherence edge weight distribution')
ylabel('Frequency')
legend('SEEG','ECoG','Location','NorthEast')


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

% figure(3);clf;
% subplot(1,2,1)
% boxplot(sw_value_ecog(:,[1,3,5]))
% ylabel('Small-world propensity')
% set(gca,'xtick',(1:3),'xticklabel',{'Broadband CC','Beta','High-Gamma'});
% ylim([0.25 1])
% title('Small-worldness of ECoG networks')
% subplot(1,2,2)
% boxplot(sw_value_seeg(:,[1,3,5]))
% set(gca,'xtick',(1:3),'xticklabel',{'Broadband CC','Beta','High-Gamma'});
% ylabel('Small-world propensity')
% title('Small-worldness of SEEG networks')
% ylim([0.25 1])

figure(3);clf;
boxplot([[sw_value_ecog(:,1);NaN],sw_value_seeg(:,1),[sw_value_ecog(:,3);NaN],...
    sw_value_seeg(:,3)])
set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
ylabel('Small-world propensity')
ylim([0.25 1.1])
h = findobj(gca,'Tag','Box');
colors = [color2; color1; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

ranksum(sw_value_ecog(:,1),sw_value_seeg(:,1))
ranksum(sw_value_ecog(:,3),sw_value_seeg(:,3))

signrank(sw_value_ecog(:,1),sw_value_ecog(:,3))

signrank(sw_value_seeg(:,1),sw_value_seeg(:,3))
%% 3B. Modularity analysis -> compute modularity & participation coefficient

[ecog_pt_modules, ecog_pt_q_vals, ecog_pt_pc, ecog_pc_all] = compute_modularity({ecog_patients.conn});
[seeg_pt_modules, seeg_pt_q_vals, seeg_pt_pc, seeg_pc_all] = compute_modularity({seeg_patients.conn});


boxplot([[ecog_pc_all(:,1);NaN],seeg_pc_all(:,1),[ecog_pc_all(:,3);NaN],...
    seeg_pc_all(:,3)])
set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
h = findobj(gca,'Tag','Box');
colors = [color2; color1; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])

ranksum(ecog_pc_all(:,1),seeg_pc_all(:,1))
ranksum(ecog_pc_all(:,3),seeg_pc_all(:,3))

%% 3C. Comparison of modules to Yeo systems & own electrode groups

% have to get 'community' labels for electrodes and Yeo systems
[ECoG_pt_purity, ECoG_module_purity] = network_similarity(ecog_pt_modules, {ecog_patients.coords}, 'systems');
[SEEG_pt_purity, SEEG_module_purity] = network_similarity(seeg_pt_modules, {seeg_patients.coords}, 'systems');

figure(1);clf
boxplot([[ECoG_module_purity(:,1);NaN],SEEG_module_purity(:,1),[ECoG_module_purity(:,3);NaN],...
    SEEG_module_purity(:,3)])
set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
h = findobj(gca,'Tag','Box');
colors = [color2; color1; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])

ranksum(ECoG_module_purity(:,1),SEEG_module_purity(:,1))
ranksum(ECoG_module_purity(:,3),SEEG_module_purity(:,3))


signrank(ECoG_module_purity(:,1),ECoG_module_purity(:,3))
signrank(SEEG_module_purity(:,1),SEEG_module_purity(:,3))

%% 4A. Node strength of the SOZ/EZ

[ecog_str_non_res, ecog_str_res, ecog_str] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'node_strength');
[seeg_str_non_res, seeg_str_res, seeg_str] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'node_strength');

figure(4);clf;
subplot(2,1,1)
boxplot([[ecog_str_non_res(:,1);NaN],[ecog_str_res(:,1);NaN],seeg_str_non_res(:,1),seeg_str_res(:,1)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score node strength')
title('Broadband CC')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([-2.5 3])
subplot(2,1,2)
boxplot([[ecog_str_non_res(:,3);NaN],[ecog_str_res(:,3);NaN],seeg_str_non_res(:,3),seeg_str_res(:,3)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score node strength')
title('Beta coherence')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([-2.5 3])

ranksum(ecog_str_non_res(:,1),ecog_str_res(:,1))
ranksum(ecog_str_non_res(:,3),ecog_str_res(:,3))

ranksum(seeg_str_non_res(:,1),seeg_str_res(:,1))
ranksum(seeg_str_non_res(:,3),seeg_str_res(:,3))

%% 4B. Betweenness centrality of the SOZ/EZ

[ecog_btw_non_res, ecog_btw_res, ecog_btw] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'betweenness_centrality');
[seeg_btw_non_res, seeg_btw_res, seeg_btw] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'betweenness_centrality');

figure(4)
subplot(2,1,1)
boxplot([[ecog_btw_non_res(:,1);NaN],[ecog_btw_res(:,1);NaN],seeg_btw_non_res(:,1),seeg_btw_res(:,1)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score betweenness centrality')
title('Broadband CC')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])
ylim([-0.5 2])
subplot(2,1,2)
boxplot([[ecog_btw_non_res(:,3);NaN],[ecog_btw_res(:,3);NaN],seeg_btw_non_res(:,3),seeg_btw_res(:,3)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score betweenness centrality')
title('Beta coherence')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])
ylim([-0.5 2])

ranksum(ecog_btw_non_res(:,1),ecog_btw_res(:,1))
ranksum(ecog_btw_non_res(:,3),ecog_btw_res(:,3))

ranksum(seeg_btw_non_res(:,1),seeg_btw_res(:,1))
ranksum(seeg_btw_non_res(:,3),seeg_btw_res(:,3))

%% 4C. Control centrality of the SOZ/EZ

[ecog_cc_non_res, ecog_cc_res, ecog_cc] = test_EZ({ecog_patients.conn},{ecog_patients.resect},'control_centrality');
[seeg_cc_non_res, seeg_cc_res, seeg_cc] = test_EZ({seeg_patients.conn},{seeg_patients.resect},'control_centrality');

figure(4)
subplot(2,1,1)
boxplot([[ecog_cc_non_res(:,1);NaN],[ecog_cc_res(:,1);NaN],seeg_cc_non_res(:,1),seeg_cc_res(:,1)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score control centrality')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])
title('Broadband CC')
ylim([-1 4])
subplot(2,1,2)
boxplot([[ecog_cc_non_res(:,3);NaN],[ecog_cc_res(:,3);NaN],seeg_cc_non_res(:,3),seeg_cc_res(:,3)])
set(gca,'xtick',(1:4),'xticklabel',{'ECoG non res','ECoG res','SEEG non res','SEEG res'});
ylabel('Z score control centrality')
h = findobj(gca,'Tag','Box');
colors = [color2; color2; color1; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 1])
title('Beta coherence')
ylim([-1 4])

ranksum(ecog_cc_non_res(:,1),ecog_cc_res(:,1))
ranksum(ecog_cc_non_res(:,3),ecog_cc_res(:,3))

ranksum(seeg_cc_non_res(:,1),seeg_cc_res(:,1))
ranksum(seeg_cc_non_res(:,3),seeg_cc_res(:,3))

%% 5A. Localization based on top metric -> All nodes -> all three methods

% node strength localization
[targets_ecog_str_raw, dice_ecog_str_raw] = localize_EZ({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.resect}, 'node_strength');
[targets_seeg_str_raw, dice_seeg_str_raw] = localize_EZ({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.resect}, 'node_strength');

% betweenness centrality localization
[targets_ecog_btw_raw, dice_ecog_btw_raw] = localize_EZ({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.resect}, 'betweenness_centrality');
[targets_seeg_btw_raw, dice_seeg_btw_raw] = localize_EZ({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.resect}, 'betweenness_centrality');

% control centrality localization
[targets_ecog_cc_raw, dice_ecog_cc_raw] = localize_EZ({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.resect}, 'control_centrality');
[targets_seeg_cc_raw, dice_seeg_cc_raw] = localize_EZ({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.resect}, 'control_centrality');

figure(6);clf;
subplot(1,3,1)
boxplot([[dice_ecog_str_raw(:,1);NaN],dice_seeg_str_raw(:,1),[dice_ecog_str_raw(:,3);NaN],dice_seeg_str_raw(:,3),[dice_ecog_str_raw(:,5);NaN],dice_seeg_str_raw(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of node strength across frequency bands')
ylim([-1 3])
subplot(1,3,2)
boxplot([[dice_ecog_btw_raw(:,1);NaN],dice_seeg_btw_raw(:,1),[dice_ecog_btw_raw(:,3);NaN],dice_seeg_btw_raw(:,3),[dice_ecog_btw_raw(:,5);NaN],dice_seeg_btw_raw(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of betweenness centrality across frequency bands')
ylim([-1 3])
subplot(1,3,3)
boxplot([[dice_ecog_cc_raw(:,1);NaN],dice_seeg_cc_raw(:,1),[dice_ecog_cc_raw(:,3);NaN],dice_seeg_cc_raw(:,3),[dice_ecog_cc_raw(:,5);NaN],dice_seeg_cc_raw(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of control centrality across frequency bands')
ylim([-1 3])

%% 5B. Where do missed localizations go? break down by target type

%% 5C. Localization based on top metric -> No WM

% generate new data
[no_wm_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'no_WM');
[no_wm_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'no_WM');

[targets_ecog_str_no_wm, dice_ecog_str_no_wm] = localize_EZ({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_ecog.resect}, 'node_strength');
[targets_seeg_str_no_wm, dice_seeg_str_no_wm] = localize_EZ({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.resect}, 'node_strength');

[targets_ecog_btw_no_wm, dice_ecog_btw_no_wm] = localize_EZ({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_ecog.resect}, 'betweenness_centrality');
[targets_seeg_btw_no_wm, dice_seeg_btw_no_wm] = localize_EZ({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.resect}, 'betweenness_centrality');

[targets_ecog_cc_no_wm, dice_ecog_cc_no_wm] = localize_EZ({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_ecog.resect}, 'control_centrality');
[targets_seeg_cc_no_wm, dice_seeg_cc_no_wm] = localize_EZ({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.resect}, 'control_centrality');

figure(7);clf;
subplot(1,3,1)
boxplot([[dice_ecog_str_no_wm(:,1);NaN],dice_seeg_str_no_wm(:,1),[dice_ecog_str_no_wm(:,3);NaN],dice_seeg_str_no_wm(:,3),[dice_ecog_str_no_wm(:,5);NaN],dice_seeg_str_no_wm(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of node strength across frequency bands')
ylim([-1 3])
subplot(1,3,2)
boxplot([[dice_ecog_btw_no_wm(:,1);NaN],dice_seeg_btw_no_wm(:,1),[dice_ecog_btw_no_wm(:,3);NaN],dice_seeg_btw_no_wm(:,3),[dice_ecog_btw_no_wm(:,5);NaN],dice_seeg_btw_no_wm(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of betweenness centralitt across frequency bands')
ylim([-1 3])
subplot(1,3,3)
boxplot([[dice_ecog_cc_no_wm(:,1);NaN],dice_seeg_cc_no_wm(:,1),[dice_ecog_cc_no_wm(:,3);NaN],dice_seeg_cc_no_wm(:,3),[dice_ecog_cc_no_wm(:,5);NaN],dice_seeg_cc_no_wm(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of node strength across frequency bands')
ylim([-1 3])

%% 5D. Localization based on top metric -> ROI reduced

% generate new data
[min_roi_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'min_ROI');
[min_roi_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'min_ROI');

[targets_ecog_str_min_roi, dice_ecog_str_min_roi] = localize_EZ({min_roi_ecog.conn}, {min_roi_ecog.coords}, {min_roi_ecog.resect}, 'node_strength');
[targets_seeg_str_min_roi, dice_seeg_str_min_roi] = localize_EZ({min_roi_seeg.conn}, {min_roi_seeg.coords}, {min_roi_seeg.resect}, 'node_strength');

[targets_ecog_btw_min_roi, dice_ecog_btw_min_roi] = localize_EZ({min_roi_ecog.conn}, {min_roi_ecog.coords}, {min_roi_ecog.resect}, 'betweenness_centrality');
[targets_seeg_btw_min_roi, dice_seeg_btw_min_roi] = localize_EZ({min_roi_seeg.conn}, {min_roi_seeg.coords}, {min_roi_seeg.resect}, 'betweenness_centrality');

[targets_ecog_cc_min_roi, dice_ecog_cc_min_roi] = localize_EZ({min_roi_ecog.conn}, {min_roi_ecog.coords}, {min_roi_ecog.resect}, 'control_centrality');
[targets_seeg_cc_min_roi, dice_seeg_cc_min_roi] = localize_EZ({min_roi_seeg.conn}, {min_roi_seeg.coords}, {min_roi_seeg.resect}, 'control_centrality');

figure(8);clf;
subplot(1,3,1)
boxplot([[dice_ecog_str_min_roi(:,1);NaN],dice_seeg_str_min_roi(:,1),[dice_ecog_str_min_roi(:,3);NaN],dice_seeg_str_min_roi(:,3),[dice_ecog_str_min_roi(:,5);NaN],dice_seeg_str_min_roi(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of node strength across frequency bands')
ylim([-1 3])
subplot(1,3,2)
boxplot([[dice_ecog_btw_min_roi(:,1);NaN],dice_seeg_btw_min_roi(:,1),[dice_ecog_btw_min_roi(:,3);NaN],dice_seeg_btw_min_roi(:,3),[dice_ecog_btw_min_roi(:,5);NaN],dice_seeg_btw_min_roi(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of betweenness centrality across frequency bands')
ylim([-1 3])
subplot(1,3,3)
boxplot([[dice_ecog_cc_min_roi(:,1);NaN],dice_seeg_cc_min_roi(:,1),[dice_ecog_cc_min_roi(:,3);NaN],dice_seeg_cc_min_roi(:,3),[dice_ecog_cc_min_roi(:,5);NaN],dice_seeg_cc_min_roi(:,5)])
set(gca,'xtick',(1:6),'xticklabel',{'ECoG Broadband','SEEG Broadband','ECoG Beta','SEEG Beta','ECoG High-gamma','SEEG High-gamma'});
ylabel('Dice coefficient')
title('Localization quality of control centrality across frequency bands')
ylim([-1 3])

%% Comparing node strength localization

% freq 1
plot_data_1 = [[dice_ecog_str_raw(:,1);NaN], dice_seeg_str_raw(:,1), dice_seeg_str_no_wm(:,1), dice_seeg_str_min_roi(:,1)];

ranksum(dice_ecog_str_raw(:,1),dice_seeg_str_raw(:,1))
ranksum(dice_seeg_str_raw(:,1),dice_ecog_str_no_wm(:,1))
ranksum(dice_seeg_str_raw(:,1),dice_ecog_str_min_roi(:,1))

% freq 3
plot_data_3 = [[dice_ecog_str_raw(:,3);NaN], dice_seeg_str_raw(:,3), dice_seeg_str_no_wm(:,3), dice_seeg_str_min_roi(:,3)];

ranksum(dice_ecog_str_raw(:,3),dice_seeg_str_raw(:,3))
ranksum(dice_seeg_str_raw(:,3),dice_ecog_str_no_wm(:,3))
ranksum(dice_seeg_str_raw(:,3),dice_ecog_str_min_roi(:,3))

%% Comparing betweenness centrality localization

% freq 1
plot_data_1 = [[dice_ecog_btw_raw(:,1);NaN], dice_seeg_btw_raw(:,1), dice_seeg_btw_no_wm(:,1), dice_seeg_btw_min_roi(:,1)]

ranksum(dice_ecog_btw_raw(:,1),dice_seeg_btw_raw(:,1))
ranksum(dice_seeg_btw_raw(:,1),dice_ecog_btw_no_wm(:,1))
ranksum(dice_seeg_btw_raw(:,1),dice_ecog_btw_min_roi(:,1))

% freq 3
plot_data_3 = [[dice_ecog_btw_raw(:,3);NaN], dice_seeg_btw_raw(:,3), dice_seeg_btw_no_wm(:,3), dice_seeg_btw_min_roi(:,3)]

ranksum(dice_ecog_btw_raw(:,3),dice_seeg_btw_raw(:,3))
ranksum(dice_seeg_btw_raw(:,3),dice_ecog_btw_no_wm(:,3))
ranksum(dice_seeg_btw_raw(:,3),dice_ecog_btw_min_roi(:,3))

%% Comparing control centrality localization

% freq 1
plot_data_1 = [[dice_ecog_cc_raw(:,1);NaN], dice_seeg_cc_raw(:,1), dice_seeg_cc_no_wm(:,1), dice_seeg_cc_min_roi(:,1)];

ranksum(dice_ecog_cc_raw(:,1),dice_seeg_cc_raw(:,1))
ranksum(dice_seeg_cc_raw(:,1),dice_ecog_cc_no_wm(:,1))
ranksum(dice_seeg_cc_raw(:,1),dice_ecog_cc_min_roi(:,1))

% freq 3
plot_data_3 = [[dice_ecog_cc_raw(:,3);NaN], dice_seeg_cc_raw(:,3), dice_seeg_cc_no_wm(:,3), dice_seeg_cc_min_roi(:,3)];

ranksum(dice_ecog_cc_raw(:,3),dice_seeg_cc_raw(:,3))
ranksum(dice_seeg_cc_raw(:,3),dice_ecog_cc_no_wm(:,3))
ranksum(dice_seeg_cc_raw(:,3),dice_ecog_cc_min_roi(:,3))

%% probably need to do distance correction....

% -> ECoG and SEEG separately -> use non-resected / ablated regions

