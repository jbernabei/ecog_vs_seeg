clear all

load HUP133_patient_struct.mat
load HUP075_patient_struct.mat

%% part 1, anatomical analysis

% 1a, renderings, with edges
% done

% 1b, top regions of electrodes and numbers in each
% HUP075:
HUP075_unique_roi = unique(HUP075_patient_struct.elec_roi);
for s = 1:length(HUP075_unique_roi)
    repeat_roi_HUP075(s) = sum(strcmp(HUP075_unique_roi{s},HUP075_patient_struct.elec_roi));
end

[B_075,I_075] = sort(repeat_roi_HUP075,'descend');
ranked_roi_HUP075 = HUP075_unique_roi(I_075)
B_075

% HUP133:
HUP133_unique_roi = unique(HUP133_patient_struct.elec_roi);
for s = 1:length(HUP133_unique_roi)
    repeat_roi_HUP133(s) = sum(strcmp(HUP133_unique_roi{s},HUP133_patient_struct.elec_roi));
end

[B_133,I_133] = sort(repeat_roi_HUP133,'descend');
ranked_roi_HUP133 = HUP133_unique_roi(I_133)
B_133

% 1c, top region - region connectivity and numbers of each

% distribution of distances
% HUP075
a = 0;
num_elecs = size(HUP075_patient_struct.mni_coords,1);
for i = 1:num_elecs
    for j = 1:num_elecs
        a = a+1;
        elec1 = HUP075_patient_struct.mni_coords(i,:);
        elec2 = HUP075_patient_struct.mni_coords(j,:);
        HUP075_inter_elec_distance(a) = sqrt(sum((elec1-elec2).^2));
    end
end

HUP075_inter_elec_distance(HUP075_inter_elec_distance==0) = [];

% HUP133
a = 0;
num_elecs = size(HUP133_patient_struct.mni_coords,1);
for i = 1:num_elecs
    for j = 1:num_elecs
        a = a+1;
        elec1 = HUP133_patient_struct.mni_coords(i,:);
        elec2 = HUP133_patient_struct.mni_coords(j,:);
        HUP133_inter_elec_distance(a) = sqrt(sum((elec1-elec2).^2));
    end
end

HUP133_inter_elec_distance(HUP133_inter_elec_distance==0) = [];

% make the figure
figure(1);clf;
subplot(1,2,1)
histogram(HUP075_inter_elec_distance)
subplot(1,2,2)
histogram(HUP133_inter_elec_distance)

skewness(HUP075_inter_elec_distance)
skewness(HUP133_inter_elec_distance)

%% 2. connectivity analysis

% 2a, distribution of connectivity strengths across bands
% HUP075
all_connectivity = [];
for k = 1:4
    all_connectivity(k).data = HUP075_patient_struct.II.mean.freq(k).data(:);
    all_connectivity(k).data(all_connectivity(k).data==0) = [];
end

figure(3);clf
subplot(2,2,1)
histogram(all_connectivity(1).data)
axis([0.05 1, 0, 3000])
title('beta')
subplot(2,2,2)
histogram(all_connectivity(2).data)
axis([0.05 1, 0, 3000])
title('low gamma')
subplot(2,2,3)
histogram(all_connectivity(3).data)
axis([0.05 1, 0, 3000])
title('high gamma')
subplot(2,2,4)
histogram(all_connectivity(4).data)
axis([0.05 1, 0, 3000])
title('broadband')

% HUP133
all_connectivity = [];
for k = 1:4
    all_connectivity(k).data = HUP133_patient_struct.II.mean.freq(k).data(:);
    all_connectivity(k).data(all_connectivity(k).data==0) = [];
end

figure(4);clf
subplot(2,2,1)
histogram(all_connectivity(1).data)
axis([0.05 1, 0, 3000])
title('beta')
subplot(2,2,2)
histogram(all_connectivity(2).data)
axis([0.05 1, 0, 3000])
title('low gamma')
subplot(2,2,3)
histogram(all_connectivity(3).data)
axis([0.05 1, 0, 3000])
title('high gamma')
subplot(2,2,4)
histogram(all_connectivity(4).data)
axis([0.05 1, 0, 3000])
title('broadband')

%% 2b, modularity analysis, with rendering
% HUP075
patient_name = 'HUP075'
subject_data = [];
for i = 1:4
    [S,Q] = modularity_und(HUP075_patient_struct.II.mean.freq(i).data,1)
    
    new_coords = HUP075_patient_struct.mni_coords;
    
    rendering_filename = sprintf('%s_mni_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = S(ch); % bad channel is red

        subject_data(ch,5) = 1; 
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)

    
    
end

% HUP133
patient_name = 'HUP133'
subject_data = [];
for i = 1:4
    [S,Q] = modularity_und(HUP133_patient_struct.II.mean.freq(i).data,1)
    
    new_coords = HUP133_patient_struct.mni_coords;
    
    rendering_filename = sprintf('%s_mni_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = S(ch); % bad channel is red

        subject_data(ch,5) = 1; 
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)

    
    
end

%% connectivity heterogeneity for each node
% HUP075
num_elecs = length(HUP075_patient_struct.elec_labels)
for i = 1:num_elecs
    for k = 1:4
    HUP075_heterogeneity(k).data(i) = var(HUP075_patient_struct.II.mean.freq(k).data(i,:));
    end
end

% HUP133
num_elecs = length(HUP133_patient_struct.elec_labels)
for i = 1:num_elecs
    for k = 1:4
    HUP133_heterogeneity(k).data(i) = var(HUP133_patient_struct.II.mean.freq(k).data(i,:));
    end
end

% make fig
figure(4);clf
subplot(2,2,1)
histogram(HUP075_heterogeneity(1).data)
subplot(2,2,2)
histogram(HUP075_heterogeneity(2).data)
subplot(2,2,3)
histogram(HUP075_heterogeneity(3).data)
subplot(2,2,4)
histogram(HUP075_heterogeneity(4).data)

figure(5);clf
subplot(2,2,1)
histogram(HUP133_heterogeneity(1).data)
subplot(2,2,2)
histogram(HUP133_heterogeneity(2).data)
subplot(2,2,3)
histogram(HUP133_heterogeneity(3).data)
subplot(2,2,4)
histogram(HUP133_heterogeneity(4).data)


%% 2d, node strength, control centrality, other metrics

for k = 1:4
    % calculate interictal node strength
    HUP075_nodestr(k).data = strengths_und(HUP075_patient_struct.II.mean.freq(k).data);
    HUP133_nodestr(k).data = strengths_und(HUP133_patient_struct.II.mean.freq(k).data);
    
    % calculate interictal control centrality
    HUP075_cc(k).data = control_centrality(HUP075_patient_struct.II.mean.freq(k).data);
    HUP133_cc(k).data = control_centrality(HUP133_patient_struct.II.mean.freq(k).data);
    
end

% plot all heterogeneity
for k = 1:4
figure(k);clf
subplot(2,2,1)
histogram(HUP075_nodestr(k).data)
title('HUP075 node strength')
subplot(2,2,2)
histogram(HUP075_cc(k).data)
title('HUP075 control centrality')
subplot(2,2,3)
histogram(HUP133_nodestr(k).data)
title('HUP133 node strength')
subplot(2,2,4)
histogram(HUP133_cc(k).data)
title('HUP133 control centrality')
end

%%
% render for HUP075
patient_name = 'HUP075'
subject_data = [];
for i = 1:4    
    new_coords = HUP075_patient_struct.mni_coords;
    
    rendering_filename = sprintf('%s_str_cc_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = HUP075_cc(i).data(ch); % bad channel is red

        subject_data(ch,5) = HUP075_nodestr(i).data(ch);
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)

    
    
end

% HUP133
patient_name = 'HUP133'
subject_data = [];
for i = 1:4
    
    new_coords = HUP133_patient_struct.mni_coords;
    
    rendering_filename = sprintf('%s_str_cc_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = HUP133_cc(i).data(ch); % bad channel is red

        subject_data(ch,5) = HUP133_nodestr(i).data(ch);
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)

    
    
end


%% 3. making the networks more similar

% 3a, remove WM

HUP075_WM_inds = find(strcmp(HUP075_patient_struct.elec_roi,'Left Cerebral White Matter'));
HUP133_WM_inds = find(strcmp(HUP133_patient_struct.elec_roi,'White_Matter'));

HUP075_mni_coords_no_WM = HUP075_patient_struct.mni_coords;
HUP075_mni_coords_no_WM(HUP075_WM_inds,:) = [];

HUP133_mni_coords_no_WM = HUP133_patient_struct.mni_coords;
HUP133_mni_coords_no_WM(HUP133_WM_inds,:) = [];

for k = 1:4
    HUP075_no_wm_II(k).data = HUP075_patient_struct.II.mean.freq(k).data;
    HUP075_no_wm_II(k).data(HUP075_WM_inds,:) = [];
    HUP075_no_wm_II(k).data(:,HUP075_WM_inds) = [];
    
    HUP133_no_wm_II(k).data = HUP133_patient_struct.II.mean.freq(k).data;
    HUP133_no_wm_II(k).data(HUP133_WM_inds,:) = [];
    HUP133_no_wm_II(k).data(:,HUP133_WM_inds) = [];
    
    figure(k);clf
    subplot(1,2,1)
    imagesc(HUP075_no_wm_II(k).data)
    subplot(1,2,2)
    imagesc(HUP133_no_wm_II(k).data)
end

%% 3b, use highest participation coefficient node in each ROI
HUP075_unique_regions = unique(HUP075_patient_struct.elec_roi);
HUP133_unique_regions = unique(HUP133_patient_struct.elec_roi);

% loop through regions to get new centroid and adjacency matrix
for k = 1:4
    for i = 1:length(HUP075_unique_regions)
        for j = 1:length(HUP075_unique_regions)
            contained_nodes_1 = find(strcmp(HUP075_patient_struct.elec_roi,HUP075_unique_regions(i)));
            contained_nodes_2 = find(strcmp(HUP075_patient_struct.elec_roi,HUP075_unique_regions(j)));
            
            HUP075_reduced_adj(k).data(i,j) = mean(mean(HUP075_patient_struct.II.mean.freq(k).data(contained_nodes_1,contained_nodes_2)));
            HUP075_reduced_coordinates(i,:) = mean(HUP075_patient_struct.mni_coords(contained_nodes_1,:));
        end
    end
end

% loop through regions to get new centroid and adjacency matrix
for k = 1:4
    for i = 1:length(HUP133_unique_regions)
        for j = 1:length(HUP133_unique_regions)
            contained_nodes_1 = find(strcmp(HUP133_patient_struct.elec_roi,HUP133_unique_regions(i)));
            contained_nodes_2 = find(strcmp(HUP133_patient_struct.elec_roi,HUP133_unique_regions(j)));
            
            HUP133_reduced_adj(k).data(i,j) = mean(mean(HUP133_patient_struct.II.mean.freq(k).data(contained_nodes_1,contained_nodes_2)));
            HUP133_reduced_coordinates(i,:) = mean(HUP133_patient_struct.mni_coords(contained_nodes_1,:),1);
        end
    end
end

for k = 1:4
    % calculate interictal node strength
    HUP075_nodestr_red(k).data = strengths_und(HUP075_reduced_adj(k).data);
    HUP133_nodestr_red(k).data = strengths_und(HUP133_reduced_adj(k).data);
    
    % calculate interictal control centrality
    HUP075_cc_red(k).data = control_centrality(HUP075_reduced_adj(k).data);
    HUP133_cc_red(k).data = control_centrality(HUP133_reduced_adj(k).data);
    
end


% render for HUP075
patient_name = 'HUP075'
subject_data = [];
for i = 1:4    
    new_coords = HUP075_reduced_coordinates;
    
    rendering_filename = sprintf('%s_str_cc_red_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = HUP075_cc_red(i).data(ch); % bad channel is red

        subject_data(ch,5) = HUP075_nodestr_red(i).data(ch);
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)

    
    
end

% HUP133
patient_name = 'HUP133'
subject_data = [];
for i = 1:4
    
    new_coords = HUP133_reduced_coordinates;
    
    rendering_filename = sprintf('%s_str_cc_red_%d.jpg',patient_name,i);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = HUP133_cc_red(i).data(ch); % bad channel is red

        subject_data(ch,5) = HUP133_nodestr_red(i).data(ch);
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_modularity.mat',rendering_filename)
    
end

%% 4. ictal analysis

% 4a, clinical, SOZ & propagation zone color coded rendering

% 4b, node strength, control centrality, etc: specificity to SOZ