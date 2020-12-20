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
color5 = [103 55 155]/255;
color6 = [78 172 91]/255;

my_colormap = make_colormap(color1,color3,color2);

% extract information on which lobes the regions of the AAL atlas belong to
lobe_table = readtable('lobes_aal.xlsx');

% extract the indices of all patients that are ECoG/SEEG respectively
ecog_patient_indices = find([hasData_field{:}] & strcmp(implant_field,'ECoG'));
seeg_patient_indices = find([hasData_field{:}] & strcmp(implant_field,'SEEG'));

% extract data into respective structures
ecog_patients = all_patients(ecog_patient_indices);
seeg_patients = all_patients(seeg_patient_indices);

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
    {ecog_patients.laterality},{ecog_patients.target}, all_inds, all_locs,lobe_table);

[SEEG_top_regions, SEEG_plot_data, SEEG_total_elecs, SEEG_loc_matrix, SEEG_bilat_score, SEEG_ipsi_focal] = ...
    rank_anatomical_targets({seeg_patients.patientID},{seeg_patients.roi},...
    {seeg_patients.laterality},{seeg_patients.target}, all_inds, all_locs,lobe_table);

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

%% for supplement:
% also make boxplots of:
% 1: intralobar
% 2: interlobar
% 3: interhemispheric
% 4: WM

% quantify intra/interlobar/interhemispheric/wm for ECoG
for pt = 1:length(ECoG_loc_matrix)
    loc_vec = ECoG_loc_matrix(pt).data(:);
    ECoG_frac_intralobar(pt,1) = sum(loc_vec==1)./sum(loc_vec~=4);
    ECoG_frac_interlobar(pt,1) = sum(loc_vec==2)./sum(loc_vec~=4);
    ECoG_frac_interhemis(pt,1) = sum(loc_vec==3)./sum(loc_vec~=4);
    ECoG_frac_WM(pt,1) = sum(loc_vec==4)./length(loc_vec);
end

% quantify intra/interlobar/interhemispheric/wm for ECoG
for pt = 1:length(SEEG_loc_matrix)
    loc_vec = SEEG_loc_matrix(pt).data(:);
    SEEG_frac_intralobar(pt,1) = sum(loc_vec==1)./sum(loc_vec~=4);
    SEEG_frac_interlobar(pt,1) = sum(loc_vec==2)./sum(loc_vec~=4);
    SEEG_frac_interhemis(pt,1) = sum(loc_vec==3)./sum(loc_vec~=4);
    SEEG_frac_WM(pt,1) = sum(loc_vec==4)./length(loc_vec);
end

% get p values for comparing ECoG/SEEG intra/interlobar/interhemis/wm
[p1, h1, stats1] = ranksum(ECoG_frac_intralobar,SEEG_frac_intralobar)
[p2, h2, stats2]  = ranksum(ECoG_frac_interlobar,SEEG_frac_interlobar)
[p3, h3, stats3]  = ranksum(ECoG_frac_interhemis,SEEG_frac_interhemis)
[p4, h4, stats4]  = ranksum(ECoG_frac_WM,SEEG_frac_WM)

e1 = mean(ECoG_frac_intralobar);
e2 = mean(ECoG_frac_interlobar);
e3 = mean(ECoG_frac_interhemis);

s1 = mean(SEEG_frac_intralobar);
s2 = mean(SEEG_frac_interlobar);
s3 = mean(SEEG_frac_interhemis);

p5 = chi2Tests([e1,e2,e3;s1,s2,s3])

% need to put GM/WM in table

[p1, p2, p3, p4]

figure(1);clf
hold on
% need to put in x labels, legend, p val bars
boxplot([[ECoG_frac_intralobar;NaN*ones(6,1)],[SEEG_frac_intralobar],...
    [ECoG_frac_interlobar;NaN*ones(6,1)],[SEEG_frac_interlobar],...
    [ECoG_frac_interhemis;NaN*ones(6,1)],[SEEG_frac_interhemis]])
h = findobj(gca,'Tag','Box');
set(gca,'xtick',[1.5:2:5.5],'xticklabel',{'Intralobar','Interlobar','Interhemispheric'});
colors = [color2; color1; color2; color1; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
legend('SEEG','ECoG','Location','NorthEast')
ylim([-0.1,1])
hold off

% fig 2
figure(2);clf;
subplot(1,2,1)
hold on
boxplot([[ECoG_ipsi_focal';NaN*ones(6,1)],SEEG_ipsi_focal']);
h = findobj(gca,'Tag','Box');
set(gca,'xtick',[1.5:2:5.5],'xticklabel',{'Intralobar','Interlobar','Interhemispheric'});
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
legend('SEEG','ECoG','Location','NorthEast')
hold off
subplot(1,2,2)
hold on
boxplot([[ECoG_bilat_score';NaN*ones(6,1)],SEEG_bilat_score']);
h = findobj(gca,'Tag','Box');
set(gca,'xtick',[1.5:2:5.5],'xticklabel',{'Intralobar','Interlobar','Interhemispheric'});
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
legend('SEEG','ECoG','Location','NorthEast')
hold off


%% 2B. Interelectrode distance comparison

[ecog_distances_all, ecog_distances_pt] = compute_interelectrode_distances({ecog_patients.coords},{ecog_patients.roi});
[seeg_distances_all, seeg_distances_pt] = compute_interelectrode_distances({seeg_patients.coords},{seeg_patients.roi});

% calculate means
m1 = mean(ecog_distances_all)
m2 = mean(seeg_distances_all)

% calculate skewness
s1 = skewness(ecog_distances_all)
s2 = skewness(seeg_distances_all)

[p1, h1, stats1]= ranksum(ecog_distances_all,seeg_distances_all)

figure(2);
hold on
histogram(ecog_distances_all,'Normalization','probability','FaceColor',color1)
histogram(seeg_distances_all,'Normalization','probability','FaceColor',color2)
ylabel('Frequency')
xlabel('Internodal distance (mm)')
title('ECoG versus SEEG internodal distances')
legend('ECoG','SEEG','Location','NorthEast')
hold off

% for supplement
% also make plots of:
% boxplots of mean internodal distance
for i = 1:length(ecog_distances_pt)
    med_ecog_dists(i) =  nanmedian(ecog_distances_pt{i}(:));
end

for i = 1:length(seeg_distances_pt)
    med_seeg_dists(i) =  nanmedian(seeg_distances_pt{i}(:));
end

[p1, h1, stats1] = ranksum(med_ecog_dists,med_seeg_dists)
figure(3);
boxplot([[med_ecog_dists';NaN*ones(6,1)],[med_seeg_dists']])
ylabel('mean internodal distance (mm)')
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
legend('SEEG','ECoG','Location','NorthEast')
ylim([35,100])
%% 2C. Connectivity analysis

% histogram of raw values / edge variability for each node -> lets do HG
% break down each histogram / bar by region (intra ROI / inter ROI)

[ecog_conn_all, ecog_conn_pt] = analyze_connectivity({ecog_patients.conn},{ecog_patients.roi});
[seeg_conn_all, seeg_conn_pt] = analyze_connectivity({seeg_patients.conn},{seeg_patients.roi});

figure(1);clf
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

m1 = mean(seeg_conn_all(:,1))
m2 = mean(ecog_conn_all(:,1))

[p,h,stats] = ranksum(ecog_conn_all(:,1),seeg_conn_all(:,1))

m3 = mean(seeg_conn_all(:,3))
m4 = mean(ecog_conn_all(:,3))

[p,h,stats] = ranksum(ecog_conn_all(:,3),seeg_conn_all(:,3))

figure(2);clf;
hold on
histogram(seeg_conn_all(:,3),'Normalization','probability','FaceColor', color2)
histogram(ecog_conn_all(:,3),'Normalization','probability','FaceColor', color1)
xlim([0.1,0.5])
title('Beta coherence edge weight distribution')
ylabel('Frequency')
legend('SEEG','ECoG','Location','NorthEast')

% for supplement
% also make plots of:
% boxplots of mean connectivity
for i = 1:length(ecog_conn_pt)
    med_ecog_conn(i) =  nanmedian(ecog_conn_pt{i}.data(:,3));
end

for i = 1:length(seeg_conn_pt)
    med_seeg_conn(i) =  nanmedian(seeg_conn_pt{i}.data(:,3));
end

[p, h, stats]= ranksum(med_ecog_conn,med_seeg_conn)
figure(3);clf;
boxplot([[med_ecog_conn';NaN*ones(6,1)],[med_seeg_conn']])
ylim([0.1 0.2])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
% 

%% do distance - connectivity regression + correction
[no_wm_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'no_WM');
[no_wm_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'no_WM');

[curve_ECoG, all_dist_ecog, all_conn_ecog,conn_data_ecog] = compute_distance_regression({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_ecog.resect});
[curve_SEEG, all_dist_seeg, all_conn_seeg,conn_data_seeg] = compute_distance_regression({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.resect});

x_axis = [1:0.1:100];
y_ecog = (curve_ECoG(3).data.a.*x_axis.^curve_ECoG(3).data.b)+curve_ECoG(3).data.c;
y_seeg = (curve_SEEG(3).data.a.*x_axis.^curve_SEEG(3).data.b)+curve_SEEG(3).data.c;

figure(1);clf
hold on
plot(all_dist_ecog(1:100:end), all_conn_ecog(1:100:end,3),'ko')
plot(all_dist_seeg(1:100:end), all_conn_seeg(1:100:end,3),'ro')
hold off

figure(2);clf;
hold on
plot(x_axis,y_ecog,'r-')
plot(x_axis,y_seeg,'b-')
hold off

%% 3A. Small worldness (or other global metric)

% !! do this without WM !!
[no_wm_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'no_WM');
[no_wm_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'no_WM');


a = 0;
for pt = 1:length(no_wm_ecog)
    a = a+1;
    for f = 1:5
        pt_adj = no_wm_ecog(pt).conn(f).data;
        [SWP_ecog,delta_C_ecog,delta_L_ecog] = small_world_propensity(pt_adj);
        sw_value_ecog(a,f) = SWP_ecog;
    end
end

a = 0;
for pt = 1:length(no_wm_seeg)
    a = a+1;
    for f = 1:5
        pt_adj = no_wm_seeg(pt).conn(f).data;
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
boxplot([[sw_value_ecog(:,3);NaN*ones(6,1)],sw_value_seeg(:,3)])
ylabel('Small-world propensity')
ylim([0.25 0.9])
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

[p,h,stats] = ranksum(sw_value_ecog(:,3),sw_value_seeg(:,3))
%% 3B. Modularity analysis -> compute modularity & participation coefficient

% !! do this without WM !!

g_val = 1;

[ecog_pt_modules, ecog_pt_q_vals, ecog_pt_pc, ecog_pc_all] = compute_modularity({no_wm_ecog.conn},g_val);
[seeg_pt_modules, seeg_pt_q_vals, seeg_pt_pc, seeg_pc_all] = compute_modularity({no_wm_seeg.conn},g_val);

figure(1);clf
boxplot([[ecog_pc_all(:,3);NaN*ones(6,1)],...
    seeg_pc_all(:,3)])
%set(gca,'xtick',(1.5:2:3.5),'xticklabel',{'Broadband CC','Beta'});
h = findobj(gca,'Tag','Box');
colors = [color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0.3 1])

mean(ecog_pc_all(:,3))

mean(seeg_pc_all(:,3))

[p,h,stats] = ranksum(ecog_pc_all(:,3),seeg_pc_all(:,3))

%% replace rest of it with distinguishability statistic
[ecog_str_D_rs] = compute_distinguishability({ecog_patients.conn},{ecog_patients.resect},'node_strength',0, {ecog_patients.roi});
[seeg_str_D_rs] = compute_distinguishability({seeg_patients.conn},{seeg_patients.resect},'node_strength',0, {seeg_patients.roi});

%%
[ecog_str_D_rs_no_WM] = compute_distinguishability({ecog_patients.conn},{ecog_patients.resect},'node_strength',1, {ecog_patients.roi});
[seeg_str_D_rs_no_WM] = compute_distinguishability({seeg_patients.conn},{seeg_patients.resect},'node_strength',1, {seeg_patients.roi});

%% reduce to min ROI
[min_roi_ecog] = modify_networks({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect}, 'min_ROI');
[min_roi_seeg] = modify_networks({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect}, 'min_ROI');

[ecog_str_D_rs_min_ROI] = compute_distinguishability({min_roi_ecog.conn}, {min_roi_ecog.resect},'node_strength',0, {ecog_patients.roi});
[seeg_str_D_rs_min_ROI] = compute_distinguishability({min_roi_seeg.conn}, {min_roi_seeg.resect},'node_strength',0, {seeg_patients.roi});


%%
freq = 3;

figure(1);clf;
subplot(1,2,1)
title('node strength')
boxplot([[ecog_str_D_rs(:,freq);NaN*ones(6,1)],seeg_str_D_rs(:,freq),...
    [ecog_str_D_rs_no_WM(:,freq);NaN*ones(6,1)],seeg_str_D_rs_no_WM(:,freq),...
    [ecog_str_D_rs_min_ROI(:,freq);NaN*ones(6,1)],seeg_str_D_rs_min_ROI(:,freq)])


%% regress out distance from connectivity matrices and transform from 0-1
% do a single regression across all patients

% ~ do we need to do all patients separately...? ~

[curve_ECoG, all_dist_ecog, all_conn_ecog,conn_data_ecog] = compute_distance_regression({no_wm_ecog.conn}, {no_wm_ecog.coords}, {no_wm_ecog.resect});
[curve_SEEG, all_dist_seeg, all_conn_seeg,conn_data_seeg] = compute_distance_regression({no_wm_seeg.conn}, {no_wm_seeg.coords}, {no_wm_seeg.resect});

[ecog_str_D_rs_distreg] = compute_distinguishability({conn_data_ecog.conn}, {no_wm_ecog.resect},'node_strength',0, {no_wm_ecog.roi});
[seeg_str_D_rs_distreg] = compute_distinguishability({conn_data_seeg.conn}, {no_wm_seeg.resect},'node_strength',0, {no_wm_seeg.roi});

%p7 = ranksum(ecog_str_D_rs_distreg(:,3),seeg_str_D_rs_distreg(:,3))
% p9 = ranksum(ecog_pc_D_rs_distreg(:,3),seeg_pc_D_rs_distreg(:,3))

%% distance correction plot
freq = 3
% figure(1);clf;
% boxplot([[ecog_str_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_str_D_rs_distreg(:,freq),...
%     [ecog_pc_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_pc_D_rs_distreg(:,freq)])

figure(1);clf;
title('node strength')
boxplot([[ecog_str_D_rs(:,freq);NaN*ones(6,1)],seeg_str_D_rs(:,freq),...
    [ecog_str_D_rs_no_WM(:,freq);NaN*ones(6,1)],seeg_str_D_rs_no_WM(:,freq),...
    [ecog_str_D_rs_min_ROI(:,freq);NaN*ones(6,1)],seeg_str_D_rs_min_ROI(:,freq),...
    [ecog_str_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_str_D_rs_distreg(:,freq)])


[p1,h,stats] = ranksum([ecog_str_D_rs(:,freq);NaN*ones(6,1)],seeg_str_D_rs(:,freq))
[p2,h,stats] = ranksum([ecog_str_D_rs_no_WM(:,freq);NaN*ones(6,1)],seeg_str_D_rs_no_WM(:,freq))
[p3,h,stats] = ranksum([ecog_str_D_rs_min_ROI(:,freq);NaN*ones(6,1)],seeg_str_D_rs_min_ROI(:,freq))
[p4,h,stats] = ranksum([ecog_str_D_rs_distreg(:,freq);NaN*ones(6,1)],seeg_str_D_rs_distreg(:,freq))

[p5,h,stats] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_no_WM(:,freq))
[p6,h,stats] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_min_ROI(:,freq))
[p7,h,stats] = signrank(ecog_str_D_rs(:,freq),ecog_str_D_rs_distreg(:,freq))

[p8,h,stats] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_no_WM(:,freq))
[p9,h,stats] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_min_ROI(:,freq))
[p10,h,stats] = signrank(seeg_str_D_rs(:,freq),seeg_str_D_rs_distreg(:,freq))

h = findobj(gca,'Tag','Box');
colors = [color2; color1;color2; color1;color2; color1;color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([-0.1 1.3])

% need to make final figure here...
% find which patients are in top third, middle third, bottom third of seeg
% in each category -> network integration, which target, surgery type
y = quantile(seeg_str_D_rs_no_WM(:,3),[0.33 0.66]);
better_inds = find(seeg_str_D_rs_no_WM(:,3)>y(2));
medium_inds = find((seeg_str_D_rs_no_WM(:,3)<y(2)).*(seeg_str_D_rs_no_WM(:,3)>y(1)));
worse_inds = find(seeg_str_D_rs_no_WM(:,3)<y(1));

a1 = seeg_pc_D_rs_no_WM(better_inds,3)
a2 = seeg_pc_D_rs_no_WM(medium_inds,3)
a3 = seeg_pc_D_rs_no_WM(worse_inds,3)

figure(2);clf;
boxplot([a1,a2,a3])

b1 = SEEG_frac_intralobar(better_inds)
b2 = SEEG_frac_intralobar(medium_inds)
b3 = SEEG_frac_intralobar(worse_inds)

c1 = SEEG_frac_interlobar(better_inds)
c2 = SEEG_frac_interlobar(medium_inds)
c3 = SEEG_frac_interlobar(worse_inds)

d1 = SEEG_frac_interhemis(better_inds)
d2 = SEEG_frac_interhemis(medium_inds)
d3 = SEEG_frac_interhemis(worse_inds)

f1 = SEEG_bilat_score(better_inds)
f2 = SEEG_bilat_score(medium_inds)
f3 = SEEG_bilat_score(worse_inds)

g1 = SEEG_ipsi_focal(better_inds)
g2 = SEEG_ipsi_focal(medium_inds)
g3 = SEEG_ipsi_focal(worse_inds)

color3 = [78 172 91]/255;
beige = [254, 249, 213]./255;
color4 = [103 55 155]/255;

figure(3);clf;
subplot(1,3,1)
boxplot([b1,b2,b3])
h = findobj(gca,'Tag','Box');
colors = [color3;beige;color4];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0.2 0.8])
subplot(1,3,2)
boxplot([c1,c2,c3])
h = findobj(gca,'Tag','Box');
colors = [color3;beige;color4];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
subplot(1,3,3)
boxplot([d1,d2,d3])
h = findobj(gca,'Tag','Box');
colors = [color3;beige;color4];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

figure(4);clf;
e1 = seeg_pc_all(better_inds,3)
e2 = seeg_pc_all(medium_inds,3)
e3 = seeg_pc_all(worse_inds,3)
boxplot([e1,e2,e3])

[p1, h1, stats1] = ranksum(b1,b2)
[p2, h2, stats2] = ranksum(b1,b3)
[p3, h3, stats3] = ranksum(b2,b3)

[p4, h1, stats1] = ranksum(d1,d2)
[p5, h2, stats2] = ranksum(d1,d3)
[p6, h3, stats3] = ranksum(d2,d3)

[p7, h1, stats1] = ranksum(f1,f2)
[p8, h2, stats2] = ranksum(f1,f3)
[p9, h3, stats3] = ranksum(f2,f3)

[p10, h1, stats1] = ranksum(g1,g2)
[p11, h2, stats2] = ranksum(g1,g3)
[p12, h3, stats3] = ranksum(g2,g3)

% do an analysis of does the EZ have high or low betweenness centrality
% compared to the rest of the network
%%
figure(5);clf;
boxplot([num_resect_s(worse_inds)',num_resect_s(medium_inds)',num_resect_s(better_inds)'])
h = findobj(gca,'Tag','Box');
colors = [color3;beige;color4];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
ylim([0 30])

[p1,h1,stats1] = ranksum(num_resect_s(worse_inds)',num_resect_s(medium_inds)')
[p2,h2,stats2] = ranksum(num_resect_s(medium_inds)',num_resect_s(better_inds)')
[p3,h3,stats3] = ranksum(num_resect_s(worse_inds)',num_resect_s(better_inds)')

%% sub-analysis by outcome
ecog_good_inds = [];
ecog_poor_inds = [];

for pt = 1:length(ecog_patients)
    if strcmp(ecog_patients(pt).outcome,'good')
        ecog_good_inds = [ecog_good_inds,pt];
    else
        ecog_poor_inds = [ecog_poor_inds,pt];
    end
end

seeg_good_inds = [];
seeg_poor_inds = [];

for pt = 1:length(seeg_patients)
    if strcmp(seeg_patients(pt).outcome,'good')
        seeg_good_inds = [seeg_good_inds,pt];
    else
        seeg_poor_inds = [seeg_poor_inds,pt];
    end
end

figure(1);clf;
subplot(1,2,1)
title('node strength')
boxplot([[ecog_str_D_rs(ecog_good_inds,freq);NaN*ones(3,1)],seeg_str_D_rs(seeg_good_inds,freq),...
    [ecog_str_D_rs_no_WM(ecog_good_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_no_WM(seeg_good_inds,freq),...
    [ecog_str_D_rs_min_ROI(ecog_good_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_min_ROI(seeg_good_inds,freq),...
    [ecog_str_D_rs_distreg(ecog_good_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_distreg(seeg_good_inds,freq)])

subplot(1,2,2)
title('node strength')
boxplot([[ecog_str_D_rs(ecog_poor_inds,freq);NaN*ones(3,1)],seeg_str_D_rs(seeg_poor_inds,freq),...
    [ecog_str_D_rs_no_WM(ecog_poor_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_no_WM(seeg_poor_inds,freq),...
    [ecog_str_D_rs_min_ROI(ecog_poor_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_min_ROI(seeg_poor_inds,freq),...
    [ecog_str_D_rs_distreg(ecog_poor_inds,freq);NaN*ones(3,1)],seeg_str_D_rs_distreg(seeg_poor_inds,freq)])


%% Hup082 adjacency matrix and 1D electrode strength and D-RS plot
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
