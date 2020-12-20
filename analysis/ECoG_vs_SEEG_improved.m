clear all

% set base path of the project
base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';

% set up workspace of the project by pulling all data and metadata
[all_patients, all_inds, all_locs, conn_field, coords_field, ...
    hasData_field, id_field, implant_field, outcome_field, resect_field, ...
    roi_field, target_field, therapy_field, region_list, region_names] = ...
    set_up_workspace_ecog_seeg(base_path);

% set up colors
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color3 = [255, 248, 209]./255; 
color4 = [103 55 155]/255;
color6 = [78 172 91]/255;

my_colormap = make_colormap(color1,color3,color2);

% extract information on which lobes the regions of the AAL atlas belong to
lobe_table = readtable('lobes_aal.xlsx');

% extract the indices of all patients that are ECoG/SEEG respectively
ecog_patient_indices = find([hasData_field{:}] & strcmp(implant_field,'ECoG'));
seeg_patient_indices = find([hasData_field{:}] & strcmp(implant_field,'SEEG'));

% extract data into respective structures
ecog_patients_raw = all_patients(ecog_patient_indices);
seeg_patients_raw = all_patients(seeg_patient_indices);

%% remove electrodes outside of brain
[ecog_patients] = remove_extraparenchymal_elecs(ecog_patients_raw);
[seeg_patients] = remove_extraparenchymal_elecs(seeg_patients_raw);

%% must plot MNI registration for everything and color nodes by resect / non-resect
% Include this in the supplement
for pt = 1:length(ecog_patients)
    adj_matrix = ecog_patients(pt).conn(1).data;
    res_bin = zeros(length(ecog_patients(pt).roi),1)
    res_bin(ecog_patients(pt).resect) = 1;
    final_elec_matrix = [ecog_patients(pt).coords,res_bin,ones(size(ecog_patients(pt).coords,1),1)];
    dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','ecog_seeg_render_2.mat',sprintf('output/render/elecs_ident_%s.jpg',ecog_patients(pt).patientID))
    delete render_elecs.node
end   

for pt = 1:length(seeg_patients)
    adj_matrix = seeg_patients(pt).conn(1).data;
    res_bin = zeros(length(seeg_patients(pt).roi),1)
    res_bin(seeg_patients(pt).resect) = 1;
    final_elec_matrix = [seeg_patients(pt).coords,res_bin,ones(size(seeg_patients(pt).coords,1),1)];
    dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','ecog_seeg_render_2.mat',sprintf('output/render/elecs_ident_%s.jpg',seeg_patients(pt).patientID))
    delete render_elecs.node
end   

%% 2A. Anatomical targets 

[ECoG_top_regions, ECoG_plot_data, ECoG_total_elecs, ECoG_loc_matrix, ECoG_bilat_score, ECoG_ipsi_focal] = ...
    rank_anatomical_targets({ecog_patients.patientID},{ecog_patients.roi},...
    {ecog_patients.laterality},{ecog_patients.target}, all_inds, all_locs, lobe_table);

[SEEG_top_regions, SEEG_plot_data, SEEG_total_elecs, SEEG_loc_matrix, SEEG_bilat_score, SEEG_ipsi_focal] = ...
    rank_anatomical_targets({seeg_patients.patientID},{seeg_patients.roi},...
    {seeg_patients.laterality},{seeg_patients.target}, all_inds, all_locs, lobe_table);

% this needs some small fixes
figure(1);clf;
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
subplot(1,2,1)
plotdata1 = ECoG_plot_data(1:15)./ECoG_total_elecs;
bh1 = bar(1:numel(plotdata1),diag(plotdata1),'stacked','FaceColor', color1);
ylabel('Fraction of total electrodes')
title('Anatomical distribution of ECoG patients');
set(gca,'xtick',(1:15),'xticklabel',ECoG_top_regions(1:15));
ylim([0, 0.11])
xtickangle(45);
subplot(1,2,2)
plotdata2 = SEEG_plot_data(1:15)./SEEG_total_elecs;
bh2 = bar(1:numel(plotdata2),diag(plotdata2),'stacked','FaceColor', color1);
bh2(6).FaceColor = color2;
bh2(9).FaceColor = color2;
bh2(13).FaceColor = color2;
bh2(15).FaceColor = color2; % contralateral coloration
title('Anatomical distribution of SEEG patients');
set(gca,'xtick',(1:15),'xticklabel',SEEG_top_regions(1:15));
ylabel('Fraction of total electrodes')
xtickangle(45);
ylim([0, 0.11])

%% 2B Bilaterality & Focality
figure(2);clf;
subplot(1,2,1)
hold on
scatter(ones(1,27),ECoG_ipsi_focal,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(ECoG_ipsi_focal) median(ECoG_ipsi_focal)],'k-','LineWidth',2)
scatter(2*ones(1,33),SEEG_ipsi_focal,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(SEEG_ipsi_focal) median(SEEG_ipsi_focal)],'k-','LineWidth',2)
xlim([0.5 2.5])
ylim([-0.2 1.2])
hold off
subplot(1,2,2)
hold on
scatter(ones(1,27),ECoG_bilat_score,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(ECoG_bilat_score) median(ECoG_bilat_score)],'k-','LineWidth',2)
scatter(2*ones(1,33),SEEG_bilat_score,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(SEEG_bilat_score) median(SEEG_bilat_score)],'k-','LineWidth',2)
xlim([0.5 2.5])
ylim([-0.2 1.2])
hold off

[p_a, h_a,stats_a] = ranksum(ECoG_ipsi_focal,SEEG_ipsi_focal)
[p_b, h_b,stats_b] = ranksum(ECoG_bilat_score,SEEG_bilat_score)

%% All network modifications:

% 1. no WM -> this becomes base that we use going forward
[no_wm_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'no_WM');
[no_wm_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'no_WM');

% 2. min ROI -> this is one point of comparison (automatically removes white matter)
[min_roi_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'min_ROI');
[min_roi_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'min_ROI');

% 3. distance regression -> corrects for effect of distance
[curve_ECoG, all_dist_ecog, all_conn_ecog,dist_reg_ecog] = compute_distance_regression({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_seeg.roi}, {no_wm_ecog.resect});
[curve_SEEG, all_dist_seeg, all_conn_seeg,dist_reg_seeg] = compute_distance_regression({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.roi}, {no_wm_seeg.resect});

%% Supplemental methods: distance - connectivity regression
x_axis1 = [1:0.1:150];
x_axis2 = [3:0.1:150];
y_ecog = (curve_ECoG(3).data.p1.*x_axis1+curve_ECoG(3).data.p2)./(x_axis1+curve_ECoG(3).data.q1);
y_seeg = (curve_SEEG(3).data.p1.*x_axis2+curve_SEEG(3).data.p2)./(x_axis2+curve_SEEG(3).data.q1);

figure(1);clf
hold on
plot(all_dist_ecog(1:200:end), all_conn_ecog(1:200:end,3),'ko')
plot(all_dist_seeg(1:200:end), all_conn_seeg(1:200:end,3),'k.')
plot(x_axis1,y_ecog,'b-','LineWidth',2)
plot(x_axis2,y_seeg,'r-','LineWidth',2)
legend('ECoG','SEEG','ECoG','SEEG')
xlabel('Internodal distance (mm)')
ylabel('Beta coherence')
hold off

%% 2B. Interelectrode distance comparison

[ecog_distances_all_base, ecog_distances_pt_base] = compute_interelectrode_distances({ecog_patients.coords},{ecog_patients.roi});
[seeg_distances_all_base, seeg_distances_pt_base] = compute_interelectrode_distances({seeg_patients.coords},{seeg_patients.roi});

[ecog_distances_all_no_WM, ecog_distances_pt_no_WM] = compute_interelectrode_distances({no_wm_ecog.coords},{no_wm_ecog.roi});
[seeg_distances_all_no_WM, seeg_distances_pt_no_WM] = compute_interelectrode_distances({no_wm_seeg.coords},{no_wm_seeg.roi});

[ecog_distances_all_min_roi, ecog_distances_pt_min_roi] = compute_interelectrode_distances({min_roi_ecog.coords},{min_roi_ecog.roi});
[seeg_distances_all_min_roi, seeg_distances_pt_min_roi] = compute_interelectrode_distances({min_roi_seeg.coords},{min_roi_seeg.roi});

[ecog_distances_all_dist_reg, ecog_distances_pt_dist_reg] = compute_interelectrode_distances({dist_reg_ecog.coords},{dist_reg_ecog.roi});
[seeg_distances_all_dist_reg, seeg_distances_pt_dist_reg] = compute_interelectrode_distances({dist_reg_seeg.coords},{dist_reg_seeg.roi});

% calculate means
m1 = mean(ecog_distances_all_no_WM);
m2 = mean(seeg_distances_all_no_WM);

m1a = mean(ecog_distances_all_base);
m2a = mean(seeg_distances_all_base);

m1b = mean(ecog_distances_all_min_roi);
m2b = mean(seeg_distances_all_min_roi);

% calculate skewness
s1 = skewness(ecog_distances_all_no_WM)
s2 = skewness(seeg_distances_all_no_WM)

s1a = skewness(ecog_distances_all_base);
s2a = skewness(seeg_distances_all_base);

s1b = skewness(ecog_distances_all_min_roi);
s2b = skewness(seeg_distances_all_min_roi);

[h1c, p1c, stats1c]= kstest2(ecog_distances_all_no_WM,seeg_distances_all_no_WM);

[h1sa, p1sa, stats1sa]= kstest2(ecog_distances_all_base,seeg_distances_all_base);
[h1sb, p1sb, stats1sb]= kstest2(ecog_distances_all_min_roi,seeg_distances_all_min_roi);

figure(1);clf
hold on
histogram(ecog_distances_all_no_WM,'Normalization','probability','FaceColor',color1)
histogram(seeg_distances_all_no_WM,'Normalization','probability','FaceColor',color2)
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG versus SEEG internodal distances')
legend('ECoG','SEEG','Location','NorthEast')
hold off

% for supplement
% also make plots of:
% boxplots of mean internodal distance
for i = 1:length(ecog_distances_pt_no_WM)
    med_ecog_dists_no_WM(i) =  nanmedian(ecog_distances_pt_no_WM{i}(:));
    med_ecog_dists_base(i) =  nanmedian(ecog_distances_pt_base{i}(:));
    med_ecog_dists_min_roi(i) =  nanmedian(ecog_distances_pt_min_roi{i}(:));
end

for i = 1:length(seeg_distances_pt_no_WM)
    med_seeg_dists_no_WM(i) =  nanmedian(seeg_distances_pt_no_WM{i}(:));
    med_seeg_dists_base(i) =  nanmedian(seeg_distances_pt_base{i}(:));
    med_seeg_dists_min_roi(i) =  nanmedian(seeg_distances_pt_min_roi{i}(:));
end

[p2, h1, stats1] = ranksum(med_ecog_dists_no_WM,med_seeg_dists_no_WM)
figure(2);clf
hold on
scatter(ones(1,27),med_ecog_dists_no_WM,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(med_ecog_dists_no_WM) median(med_ecog_dists_no_WM)],'k-','LineWidth',2)
scatter(2*ones(1,33),med_seeg_dists_no_WM,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(med_seeg_dists_no_WM) median(med_seeg_dists_no_WM)],'k-','LineWidth',2)
xlim([0.5 2.5])
ylim([35 100])
hold off

%% Extra figures: visualize and stats for other bands
figure(3);clf
subplot(1,2,1)
hold on
histogram(ecog_distances_all_base,'Normalization','probability','FaceColor',color1)
histogram(seeg_distances_all_base,'Normalization','probability','FaceColor',color2)
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG versus SEEG internodal distances')
legend('ECoG','SEEG','Location','NorthEast')
hold off
subplot(1,2,2)
hold on
histogram(ecog_distances_all_min_roi,'Normalization','probability','FaceColor',color1)
histogram(seeg_distances_all_min_roi,'Normalization','probability','FaceColor',color2)
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG versus SEEG internodal distances')
legend('ECoG','SEEG','Location','NorthEast')
hold off

% now patient level
figure(4);clf
[p2a, h1, stats1] = ranksum(med_ecog_dists_base,med_seeg_dists_base);
[p2b, h1, stats1] = ranksum(med_ecog_dists_min_roi,med_seeg_dists_min_roi);
subplot(1,2,1)
hold on
boxplot([[med_ecog_dists_base';NaN*ones(6,1)],[med_seeg_dists_base']])
ylabel('mean internodal distance (mm)')
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
hold off
subplot(1,2,2)
hold on
boxplot([[med_ecog_dists_min_roi';NaN*ones(6,1)],[med_seeg_dists_min_roi']])
ylabel('mean internodal distance (mm)')
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
hold off
p2a
p2b

%% 3A. Connectivity analysis

% histogram of raw values / edge variability for each node -> lets do HG
% break down each histogram / bar by region (intra ROI / inter ROI)

[ecog_conn_all, ecog_conn_pt] = analyze_connectivity({ecog_patients.conn},{ecog_patients.roi});
[seeg_conn_all, seeg_conn_pt] = analyze_connectivity({seeg_patients.conn},{seeg_patients.roi});

[ecog_conn_all_no_WM, ecog_conn_pt_no_WM] = analyze_connectivity({no_wm_ecog.conn},{no_wm_ecog.roi});
[seeg_conn_all_no_WM, seeg_conn_pt_no_WM] = analyze_connectivity({no_wm_seeg.conn},{no_wm_seeg.roi});

[ecog_conn_all_min_roi, ecog_conn_pt_min_roi] = analyze_connectivity({min_roi_ecog.conn},{min_roi_ecog.roi});
[seeg_conn_all_min_roi, seeg_conn_pt_min_roi] = analyze_connectivity({min_roi_seeg.conn},{min_roi_seeg.roi});

figure(1);clf
hold on
histogram(seeg_conn_all_no_WM(:,3),'Normalization','probability','FaceColor', color2)
histogram(ecog_conn_all_no_WM(:,3),'Normalization','probability','FaceColor', color1)
xlim([0.1,0.6])
title('Beta coherence edge weight distribution')
ylabel('Frequency')
legend('SEEG','ECoG','Location','NorthEast')

m3 = mean(seeg_conn_all_no_WM(:,3))
m4 = mean(ecog_conn_all_no_WM(:,3))

[p,h,stats] = ranksum(ecog_conn_all_no_WM(:,3),seeg_conn_all_no_WM(:,3))

% also make plots of:
% boxplots of mean connectivity
for i = 1:length(ecog_conn_pt_no_WM)
    med_ecog_conn_no_WM(i) =  nanmedian(ecog_conn_pt_no_WM{i}.data(:,3));
    med_ecog_conn_base(i) =  nanmedian(ecog_conn_pt{i}.data(:,3));
    med_ecog_conn_min_roi(i) =  nanmedian(ecog_conn_pt_min_roi{i}.data(:,3));
end

for i = 1:length(seeg_conn_pt)
    med_seeg_conn_no_WM(i) =  nanmedian(seeg_conn_pt_no_WM{i}.data(:,3));
    med_seeg_conn_base(i) =  nanmedian(seeg_conn_pt{i}.data(:,3));
    med_seeg_conn_min_roi(i) =  nanmedian(seeg_conn_pt_min_roi{i}.data(:,3));
end

[p3, h, stats]= ranksum(med_ecog_conn_no_WM,med_seeg_conn_no_WM)
figure(3);clf;
hold on
scatter(ones(1,27),med_ecog_conn_no_WM,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(med_ecog_conn_no_WM) median(med_ecog_conn_no_WM)],'k-','LineWidth',2)
scatter(2*ones(1,33),med_seeg_conn_no_WM,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(med_seeg_conn_no_WM) median(med_seeg_conn_no_WM)],'k-','LineWidth',2)
xlim([0.5 2.5])
ylim([0.1 0.2])

%% additional for supplement
figure(4);clf
subplot(1,2,1)
[p3a, h, stats]= ranksum(med_ecog_conn_base,med_seeg_conn_base);
boxplot([[med_ecog_conn_base';NaN*ones(6,1)],[med_seeg_conn_base']])
ylim([0.1 0.2])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

subplot(1,2,2)
[p3b, h, stats]= ranksum(med_ecog_conn_min_roi,med_seeg_conn_min_roi);
boxplot([[med_ecog_conn_min_roi';NaN*ones(6,1)],[med_seeg_conn_min_roi']])
ylim([0.1 0.2])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
p3a
p3b

%% 3B. small-worldness

a = 0;
for pt = 1:length(no_wm_ecog)
    a = a+1;
    for f = 3
        % no WM
        pt_adj = no_wm_ecog(pt).conn(f).data;
        [SWP_ecog,delta_C_ecog,delta_L_ecog] = small_world_propensity(pt_adj);
        sw_value_ecog(a,f) = SWP_ecog;
        
        % base
        pt_adj = ecog_patients(pt).conn(f).data;
        [SWP_ecog,delta_C_ecog,delta_L_ecog] = small_world_propensity(pt_adj);
        sw_value_ecog_base(a,f) = SWP_ecog;
        
        % min-roi
        pt_adj = min_roi_ecog(pt).conn(f).data;
        [SWP_ecog,delta_C_ecog,delta_L_ecog] = small_world_propensity(pt_adj);
        sw_value_ecog_min_roi(a,f) = SWP_ecog;
        
    end
end

a = 0;
for pt = 1:length(no_wm_seeg)
    a = a+1;
    for f = 3
        % no WM
        pt_adj = no_wm_seeg(pt).conn(f).data;
        [SWP_seeg,delta_C_seeg,delta_L_seeg] = small_world_propensity(pt_adj);
        sw_value_seeg(a,f) = SWP_seeg;
        
        % base
        pt_adj = seeg_patients(pt).conn(f).data;
        [SWP_seeg,delta_C_seeg,delta_L_seeg] = small_world_propensity(pt_adj);
        sw_value_seeg_base(a,f) = SWP_seeg;
        
        % min-roi
        pt_adj = min_roi_seeg(pt).conn(f).data;
        [SWP_seeg,delta_C_seeg,delta_L_seeg] = small_world_propensity(pt_adj);
        sw_value_seeg_min_roi(a,f) = SWP_seeg;
    end
end

figure(3);clf;
hold on
scatter(ones(1,27),sw_value_ecog(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(sw_value_ecog(:,3)) median(sw_value_ecog(:,3))],'k-','LineWidth',2)
scatter(2*ones(1,33),sw_value_seeg(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(sw_value_seeg(:,3)) median(sw_value_seeg(:,3))],'k-','LineWidth',2)
xlim([0.5 2.5])
ylabel('Small-world propensity')
ylim([0.25 0.9])
hold off

[p_4,h_4,stats_4] = ranksum(sw_value_ecog(:,3),sw_value_seeg(:,3))

%% for supplement
figure(4);clf;
subplot(1,2,1)
boxplot([[sw_value_ecog_base(:,3);NaN*ones(6,1)],sw_value_seeg_base(:,3)])
ylabel('Small-world propensity')
ylim([0.25 0.9])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

[p4a,h,stats] = ranksum(sw_value_ecog_base(:,3),sw_value_seeg_base(:,3));
p4a

subplot(1,2,2)
boxplot([[sw_value_ecog_min_roi(:,3);NaN*ones(6,1)],sw_value_seeg_min_roi(:,3)])
ylabel('Small-world propensity')
ylim([0.25 0.9])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

[p4b,h,stats] = ranksum(sw_value_ecog_min_roi(:,3),sw_value_seeg_min_roi(:,3));
p4b
%% 3B. Modularity analysis -> compute modularity & participation coefficient

% !! do this without WM !!

g_val = 1;

[ecog_pt_modules, ecog_pt_q_vals, ecog_pt_pc, ecog_pc_all_no_WM] = compute_modularity({no_wm_ecog.conn},g_val);
[seeg_pt_modules, seeg_pt_q_vals, seeg_pt_pc, seeg_pc_all_no_WM] = compute_modularity({no_wm_seeg.conn},g_val);

[ecog_pt_modules_base, ecog_pt_q_vals, ecog_pt_pc_base, ecog_pc_all_base] = compute_modularity({ecog_patients.conn},g_val);
[seeg_pt_modules_base, seeg_pt_q_vals, seeg_pt_pc_base, seeg_pc_all_base] = compute_modularity({seeg_patients.conn},g_val);

[ecog_pt_modules_min_roi, ecog_pt_q_vals, ecog_pt_pc_min_roi, ecog_pc_all_min_roi] = compute_modularity({min_roi_ecog.conn},g_val);
[seeg_pt_modules_min_roi, seeg_pt_q_vals, seeg_pt_pc_min_roi, seeg_pc_all_min_roi] = compute_modularity({min_roi_seeg.conn},g_val);


figure(1);clf
hold on
scatter(ones(1,27),ecog_pc_all_no_WM(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(ecog_pc_all_no_WM(:,3)) median(ecog_pc_all_no_WM(:,3))],'k-','LineWidth',2)
scatter(2*ones(1,33),seeg_pc_all_no_WM(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(seeg_pc_all_no_WM(:,3)) median(seeg_pc_all_no_WM(:,3))],'k-','LineWidth',2)
xlim([0.5 2.5])
ylim([0.3 1])
hold off

mean(ecog_pc_all_no_WM(:,3));
mean(seeg_pc_all_no_WM(:,3));
[p5,h,stats] = ranksum(ecog_pc_all_no_WM(:,3),seeg_pc_all_no_WM(:,3));
p5

%% for supplement
figure(2);clf;
subplot(1,2,1)
boxplot([[ecog_pc_all_base(:,3);NaN*ones(6,1)],...
    seeg_pc_all_base(:,3)])
%set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

mean(ecog_pc_all_base(:,3));
mean(seeg_pc_all_base(:,3));
[p5a,h,stats] = ranksum(ecog_pc_all_base(:,3),seeg_pc_all_base(:,3));
p5a

subplot(1,2,2)
boxplot([[ecog_pc_all_min_roi(:,3);NaN*ones(6,1)],...
    seeg_pc_all_min_roi(:,3)])
%set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

mean(ecog_pc_all_min_roi(:,3));
mean(seeg_pc_all_min_roi(:,3));
[p5b,h,stats] = ranksum(ecog_pc_all_min_roi(:,3),seeg_pc_all_min_roi(:,3));
p5b

%% compute distinguishability statistic under different conditions

[ecog_str_D_rs] = compute_distinguishability({ecog_patients.conn},{ecog_patients.resect},'node_strength',0, {ecog_patients.roi});
[seeg_str_D_rs] = compute_distinguishability({seeg_patients.conn},{seeg_patients.resect},'node_strength',0, {seeg_patients.roi});

[ecog_str_D_rs_no_WM] = compute_distinguishability({ecog_patients.conn},{ecog_patients.resect},'node_strength',1, {ecog_patients.roi});
[seeg_str_D_rs_no_WM] = compute_distinguishability({seeg_patients.conn},{seeg_patients.resect},'node_strength',1, {seeg_patients.roi});

[ecog_str_D_rs_min_ROI] = compute_distinguishability({min_roi_ecog.conn}, {min_roi_ecog.resect},'node_strength',0, {ecog_patients.roi});
[seeg_str_D_rs_min_ROI] = compute_distinguishability({min_roi_seeg.conn}, {min_roi_seeg.resect},'node_strength',0, {seeg_patients.roi});

[ecog_str_D_rs_distreg] = compute_distinguishability({dist_reg_ecog.conn}, {no_wm_ecog.resect},'node_strength',0, {no_wm_ecog.roi});
[seeg_str_D_rs_distreg] = compute_distinguishability({dist_reg_seeg.conn}, {no_wm_seeg.resect},'node_strength',0, {no_wm_seeg.roi});

%% distance correction plot
freq = 3;
figure(1);clf;
title('node strength')
% boxplot([[ecog_str_D_rs(:,freq);NaN*ones(6,1)],seeg_str_D_rs(:,freq),...
%     [ecog_str_D_rs_no_WM(:,freq);NaN*ones(6,1)],seeg_str_D_rs_no_WM(:,freq),...
%     [ecog_str_D_rs_min_ROI(:,freq);NaN*ones(6,1)],seeg_str_D_rs_min_ROI(:,freq),...
%     [ecog_str_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_str_D_rs_distreg(:,freq)])
hold on
scatter(ones(1,27),ecog_str_D_rs(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(ecog_str_D_rs(:,3)) median(ecog_str_D_rs(:,3))],'k-','LineWidth',2)
scatter(2*ones(1,33),seeg_str_D_rs(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(seeg_str_D_rs(:,3)) median(seeg_str_D_rs(:,3))],'k-','LineWidth',2)
scatter(3*ones(1,27),ecog_str_D_rs_no_WM(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([2.75 3.25],[median(ecog_str_D_rs_no_WM(:,3)) median(ecog_str_D_rs_no_WM(:,3))],'k-','LineWidth',2)
scatter(4*ones(1,33),seeg_str_D_rs_no_WM(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([3.75 4.25],[median(seeg_str_D_rs_no_WM(:,3)) median(seeg_str_D_rs_no_WM(:,3))],'k-','LineWidth',2)
scatter(5*ones(1,27),ecog_str_D_rs_min_ROI(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([4.75 5.25],[median(ecog_str_D_rs_min_ROI(:,3)) median(ecog_str_D_rs_min_ROI(:,3))],'k-','LineWidth',2)
scatter(6*ones(1,33),seeg_str_D_rs_min_ROI(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([5.75 6.25],[median(seeg_str_D_rs_min_ROI(:,3)) median(seeg_str_D_rs_min_ROI(:,3))],'k-','LineWidth',2)
scatter(7*ones(1,27),ecog_str_D_rs_distreg(:,3),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([6.75 7.25],[median(ecog_str_D_rs_distreg(:,3)) median(ecog_str_D_rs_distreg(:,3))],'k-','LineWidth',2)
scatter(8*ones(1,33),seeg_str_D_rs_distreg(:,3),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([7.75 8.25],[median(seeg_str_D_rs_distreg(:,3)) median(seeg_str_D_rs_distreg(:,3))],'k-','LineWidth',2)
hold off

[p1,h,stats1] = ranksum([ecog_str_D_rs(:,freq);NaN*ones(6,1)],seeg_str_D_rs(:,freq));
[p2,h,stats2] = ranksum([ecog_str_D_rs_no_WM(:,freq);NaN*ones(6,1)],seeg_str_D_rs_no_WM(:,freq));
[p3,h,stats3] = ranksum([ecog_str_D_rs_min_ROI(:,freq);NaN*ones(6,1)],seeg_str_D_rs_min_ROI(:,freq));
[p4,h,stats4] = ranksum([ecog_str_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_str_D_rs_distreg(:,freq));

[p5,h,stats5] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_no_WM(:,freq));
[p6,h,stats6] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_min_ROI(:,freq));
[p7,h,stats7] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_distreg(:,freq));

[p8,h,stats8] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_no_WM(:,freq));
[p9,h,stats9] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_min_ROI(:,freq));
[p10,h,stats10] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_distreg(:,freq));


ylim([-0.1 1.3])
[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10]

%% need to make final figure here. Why does model perform well in some SEEG patients and not others?
seeg_high_drs = find((seeg_str_D_rs_no_WM(:,3) > 0.66) +  (seeg_str_D_rs_no_WM(:,3) < 0.33));
seeg_low_drs = find((seeg_str_D_rs_no_WM(:,3) < 0.66).*(seeg_str_D_rs_no_WM(:,3) > 0.33));

% eight potential tests:
% bilaterality
bilat_high = SEEG_bilat_score(seeg_high_drs);
bilat_low = SEEG_bilat_score(seeg_low_drs);
% figure(1);clf;
% hold on
% scatter(ones(1,22),bilat_low,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
% plot([0.75 1.25],[median(bilat_low) median(bilat_low)],'k-','LineWidth',2)
% scatter(2*ones(1,11),bilat_high,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
% plot([1.75 2.25],[median(bilat_high) median(bilat_high)],'k-','LineWidth',2)
% hold off
p1 = ranksum(bilat_high,bilat_low)

% unilateral focal
ipsi_high = SEEG_ipsi_focal(seeg_high_drs);
ipsi_low = SEEG_ipsi_focal(seeg_low_drs);

p2 = ranksum(ipsi_high,ipsi_low)

% outcome
seeg_outcome = {seeg_patients.outcome}';
high_good = sum(strcmp(seeg_outcome(seeg_high_drs),'good'));
high_poor = sum(strcmp(seeg_outcome(seeg_high_drs),'poor'));
low_good = sum(strcmp(seeg_outcome(seeg_low_drs),'good'));
low_poor = sum(strcmp(seeg_outcome(seeg_low_drs),'poor'));

p3 = chi2Tests([high_good,high_poor;low_good,low_poor])

% # nodes targeted/therapy
resect_high = num_resect_s(seeg_high_drs);
resect_low = num_resect_s(seeg_low_drs);

[p4, h4, stats4] = ranksum(resect_high,resect_low)
% figure(1);clf;
% hold on
% scatter(ones(1,27),resect_low,'MarkerEdgeColor',color3,'MarkerFaceColor',color3,'jitter','on')
% plot([0.75 1.25],[median(resect_low) median(resect_low)],'k-','LineWidth',2)
% scatter(2*ones(1,6),resect_high,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
% plot([1.75 2.25],[median(resect_high) median(resect_high)],'k-','LineWidth',2)
% ylim([0 30])
% hold off

%primary target
seeg_target = {seeg_patients.target}';
high_temp = sum(strcmp(seeg_target(seeg_high_drs),'Temporal'));
high_front = sum(strcmp(seeg_target(seeg_high_drs),'Frontal'));
low_temp = sum(strcmp(seeg_target(seeg_low_drs),'Temporal'));
low_front = sum(strcmp(seeg_target(seeg_low_drs),'Frontal'));
high_par = sum(strcmp(seeg_target(seeg_high_drs),'Parietal'));
high_ins = sum(strcmp(seeg_target(seeg_high_drs),'Insular'));
low_par = sum(strcmp(seeg_target(seeg_low_drs),'Parietal'));
low_ins = sum(strcmp(seeg_target(seeg_low_drs),'Insular'));

p5 = chi2Tests([high_temp,high_front,high_par,high_ins;low_temp,low_front,low_par,low_ins])

% total nodes
nodes_high = num_nodes_s(seeg_high_drs);
nodes_low = num_nodes_s(seeg_low_drs);

[p6, h6, stats6] = ranksum(nodes_high,nodes_low)

% network integration
pc_high = seeg_pc_all_no_WM(seeg_high_drs,3);
pc_low = seeg_pc_all_no_WM(seeg_low_drs,3);

figure(1);clf;
hold on
scatter(ones(1,17)',pc_low,'MarkerEdgeColor',color6,'MarkerFaceColor',color6,'jitter','on')
plot([0.75 1.25],[median(pc_low) median(pc_low)],'k-','LineWidth',2)
scatter(2*ones(1,16)',pc_high,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
plot([1.75 2.25],[median(pc_high) median(pc_high)],'k-','LineWidth',2)
hold off

p7 = ranksum(pc_high,pc_low)

% small worldness
sw_high = sw_value_seeg(seeg_high_drs,3);
sw_low = sw_value_seeg(seeg_low_drs,3);

p8 = ranksum(sw_high,sw_low)

%therapy
seeg_therapy = {seeg_patients.therapy}';
high_abl = sum(strcmp(seeg_therapy(seeg_high_drs),'Ablation'));
high_res = sum(strcmp(seeg_therapy(seeg_high_drs),'Resection'));
low_abl = sum(strcmp(seeg_therapy(seeg_low_drs),'Ablation'));
low_res = sum(strcmp(seeg_therapy(seeg_low_drs),'Resection'));

p9 = chi2Tests([high_abl,high_res;low_abl,low_res])

%%
figure(1);clf;imagesc(corr([seeg_str_D_rs(:,freq),seeg_str_D_rs_no_WM(:,freq),seeg_str_D_rs_min_ROI(:,freq),seeg_str_D_rs_distreg(:,freq)]))
colormap(color_bar)
colorbar
%% metadata for table
for pt = 1:length(ecog_patients)
    num_nodes(pt) = length(ecog_patients(pt).roi);
    num_WM(pt) = sum(ecog_patients(pt).roi==9171);
    num_resect(pt) = length(ecog_patients(pt).resect);
end

for pt = 1:length(seeg_patients)
    num_nodes_s(pt) = length(seeg_patients(pt).roi);
    num_WM_s(pt) = sum(seeg_patients(pt).roi==9171);
    num_resect_s(pt) = length(seeg_patients(pt).resect);
end

[mean(num_nodes), std(num_nodes), mean(num_nodes_s), std(num_nodes_s), ranksum(num_nodes,num_nodes_s);
 mean(num_WM),    std(num_WM),    mean(num_WM_s),    std(num_WM_s),    ranksum(num_WM,num_WM_s);
 mean(num_resect),std(num_resect),mean(num_resect_s),std(num_resect_s),ranksum(num_resect,num_resect_s);]

%% make graphics
load color_bar

figure(1);clf
imagesc(ecog_patients(10).conn(3).data)
colormap(color_bar)
colorbar

figure(2);clf
imagesc(sum(ecog_patients(10).conn(3).data,2))
colormap(color_bar)

% make scatter for resected
res_bin = zeros(length(ecog_patients(10).roi),1);
res_bin(ecog_patients(10).resect) = 1;

nodestr = sum(ecog_patients(10).conn(3).data);

figure(3);clf;
hold on
scatter(ones(1,length(nodestr(res_bin==0))),nodestr(res_bin==0),'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
scatter(1.5*ones(1,length(nodestr(res_bin==1))),nodestr(res_bin==1),'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
xlim([0.5 2])

%% make renderings of same (seeg) patient under regular, no WM, and min roi
% conditions for supplemental figure try HUP133 6

% full matrix
adj_matrix = seeg_patients(6).conn(3).data;
res_bin = zeros(length(seeg_patients(6).roi),1)
res_bin(seeg_patients(6).resect) = 1;
save('test.edge','adj_matrix','-ascii');
final_elec_matrix = [seeg_patients(6).coords,res_bin,ones(size(seeg_patients(6).coords,1),1)];
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','test.edge','ecog_seeg_render_3.mat','output/figure/HUP133_all_nodes.jpg')
delete render_elecs.node 
delete test.edge 

% no WM
adj_matrix = no_wm_seeg(6).conn(3).data;
res_bin = zeros(length(no_wm_seeg(6).roi),1)
res_bin(no_wm_seeg(6).resect) = 1;
save('test.edge','adj_matrix','-ascii');
final_elec_matrix = [no_wm_seeg(6).coords,res_bin,ones(size(no_wm_seeg(6).coords,1),1)];
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','test.edge','ecog_seeg_render_3.mat','output/figure/HUP133_no_wm.jpg')
delete render_elecs.node 
delete test.edge 

% min ROI
adj_matrix = min_roi_seeg(6).conn(3).data;
res_bin = zeros(length(min_roi_seeg(6).roi),1)
res_bin(min_roi_seeg(6).resect) = 1;
save('test.edge','adj_matrix','-ascii');
final_elec_matrix = [min_roi_seeg(6).coords,res_bin,ones(size(min_roi_seeg(6).coords,1),1)];
dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','test.edge','ecog_seeg_render_3.mat','output/figure/HUP133_min_roi.jpg')
delete render_elecs.node 
delete test.edge 

%% 
figure(1);clf
imagesc(no_wm_seeg(10).conn(3).data)
colormap(color_bar)
colorbar

figure(2);clf
imagesc(dist_reg_seeg(10).conn(3).data)
caxis([-0.5,0.5])
colormap(color_bar)
colorbar
