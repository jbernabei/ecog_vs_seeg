%% ECoG SEEG network comparison analysis script
% John Bernabei
% email: johnbe at seas dot upenn dot edu
% MD/PhD Student, Litt Laboratory
% Center for neuroengineering & therapeutics
% Summer 2020

% The purpose of this script is to analyze data and generate figures for a
% project comparing the networks of patients with drug-resistant epilepsy
% implanted with two distinct paradigms of intracranial electrodes:
% subdural grids & strips +/- depth electrodes (ECoG) and
% stereoelectroencephalography (SEEG) in which only depth electrodes are used.

% Dependencies: 
%               BrainNet viewer: for renderings
%               Brain connectivity toolbox: for network analyses

% all other functions are contained in the package

% Input: patient data structures (see example) which contains adjacency
% matrices, electrode label, coordinate, and ROI information.

%% Set up workspace and define paired patient data

clear all

load seeg_struct
load ecog_struct
         
base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_vs_seeg/';
group_figures_path = '/Users/jbernabei/Documents/PhD_Research/ecog_vs_seeg/figures/';

  
right_left_ecog = [0 0 1 0 0 1 1 0];
right_left_seeg = [1 0 0 0 0 0 1 1];
%% Figure 1A: anatomical renderings

% Loop through all patients
for s = 1:8
    
    % Load in patient ID
    ID_ecog = ecog_struct(s).ID;
    ID_seeg = seeg_struct(s).ID;
    
    % Set output path
    output_ecog = sprintf('%s/',strcat(base_path, ID_ecog));
    output_seeg = sprintf('%s/',strcat(base_path, ID_seeg));
    
    % Set filename
    ecog_filename = sprintf('base_rendering_%s.jpg',ID_ecog);
    seeg_filename = sprintf('base_rendering_%s.jpg',ID_seeg);
    
    % Get mni coordinates
    mni_ecog = ecog_struct(s).mni_coords;
    mni_seeg = seeg_struct(s).mni_coords;
    
    % Get color data
    color_ecog = -1*ones(size(mni_ecog));
    color_seeg = -1*ones(size(mni_seeg));
    
    % Get color data
    size_ecog = ones(size(mni_ecog));
    size_seeg = ones(size(mni_seeg));
    
    % Create rendering
    create_rendering(output_ecog, ecog_filename, mni_ecog, color_ecog, size_ecog)
    create_rendering(output_seeg, seeg_filename, mni_seeg, color_seeg, size_seeg)
    
end

%% Figure 2A: distribution of inter-electrode distances
all_ecog_distance = []
all_seeg_distance = []

for s = 1:8
    % Get mni coordinates
    mni_ecog = ecog_struct(s).mni_coords;
    mni_seeg = seeg_struct(s).mni_coords;
    
    % number of electrodes in each
    num_elec_ecog = size(mni_ecog,1);
    num_elec_seeg = size(mni_seeg,1);
    
    a = 0;
    % loop through all for ecog
    for i = 1:num_elec_ecog
        for j = 1:num_elec_ecog
            a = a+1;
            elec1 = mni_ecog(i,:);
            elec2 = mni_ecog(j,:);
            ecog_distance(s).data(i,j) = sqrt(sum((elec1-elec2).^2));
            for k = 1:5
                ecog_connectivity(s).data(i,j,k) = ecog_struct(s).ii(k).data(i,j);
            end
        end
    end
    
    all_ecog_distance = [all_ecog_distance; ecog_distance(s).data(:)];
    
    a = 0;
    % loop through all for seeg
    for i = 1:num_elec_seeg
        for j = 1:num_elec_seeg
            a = a+1;
            elec1 = mni_seeg(i,:);
            elec2 = mni_seeg(j,:);
            seeg_distance(s).data(i,j) = sqrt(sum((elec1-elec2).^2));
            for k = 1:5
                seeg_connectivity(s).data(i,j,k) = seeg_struct(s).ii(k).data(i,j);
            end
        end
    end
    
    all_seeg_distance = [all_seeg_distance; seeg_distance(s).data(:)];
end

all_ecog_distance(find(all_ecog_distance==0)) = [];
all_seeg_distance(find(all_seeg_distance==0)) = [];

s1 = skewness(all_ecog_distance)
s2 = skewness(all_seeg_distance)

% make figure
figure(1);clf;
subplot(1,2,1)
histogram(all_ecog_distance)
xlabel('Internodal distance (mm)')
ylabel('Number of edges')
title('Group - level ECoG')
subplot(1,2,2)
histogram(all_seeg_distance)
xlabel('Internodal distance (mm)')
ylabel('Number of edges')
title('Group - level SEEG')

%% Find brain ROI breakdown

for s = 1:8
    % loop through brain regions ECog
        if right_left_ecog(s)==1
            ipsilateral = contains(ecog_struct(s).elec_roi,'_R');
            contralateral = contains(ecog_struct(s).elec_roi,'_L');
        else
            ipsilateral = contains(ecog_struct(s).elec_roi,'_L');
            contralateral = contains(ecog_struct(s).elec_roi,'_R');
        end
     
    %
    for r = 1:length(ecog_struct(s).elec_roi)
        processed_roi = 'n/a'
        
        roi_name = ecog_struct(s).elec_roi{r};
        
        new_roi = strip(roi_name,'right','R');
        new_roi = strip(new_roi,'right','L');
        
        if ipsilateral(r)==1
            processed_roi = strcat(new_roi,'ipsilateral');
        elseif contralateral(r)==1
            processed_roi = strcat(new_roi,'contralateral');
        end
        
        all_processed_roi(s).ecog{r} = processed_roi;
        
    end
    
    
    
    % loop through brain regions SEEG
        if right_left_seeg(s)==1
            ipsilateral = contains(seeg_struct(s).elec_roi,'_R');
            contralateral = contains(seeg_struct(s).elec_roi,'_L');
        else
            ipsilateral = contains(seeg_struct(s).elec_roi,'_L');
            contralateral = contains(seeg_struct(s).elec_roi,'_R');
        end
        
    for r = 1:length(seeg_struct(s).elec_roi)
        processed_roi = 'n/a'
        
        roi_name = seeg_struct(s).elec_roi{r};
        
        new_roi = strip(roi_name,'right','R');
        new_roi = strip(new_roi,'right','L');
        
        if ipsilateral(r)==1
            processed_roi = strcat(new_roi,'ipsilateral');
        elseif contralateral(r)==1
            processed_roi = strcat(new_roi,'contralateral');
        end
        
        all_processed_roi(s).seeg{r} = processed_roi;
    end
end

% Now we have to go through all patients and determine ranked list of ROI
% overall in ECoG and apply that list to SEEG
all_roi_ECoG = []
all_roi_SEEG = []

for s = 1:8
    all_roi_ECoG = [all_roi_ECoG; all_processed_roi(s).ecog'];
    all_roi_SEEG = [all_roi_SEEG; all_processed_roi(s).seeg'];
end

unique_regions_ECoG = unique(all_roi_ECoG);
unique_regions_SEEG = unique(all_roi_SEEG);

for r = 1:length(unique_regions_ECoG)
    freq_region_ecog(r) = sum(strcmp(all_roi_ECoG,unique_regions_ECoG{r}));
end

for r = 1:length(unique_regions_SEEG)
    freq_region_seeg(r) = sum(strcmp(all_roi_SEEG,unique_regions_SEEG{r}));
end

% Now find top 25 for ECoG
[b1, i1] = sort(freq_region_ecog,'descend');
top_regions_ecog = unique_regions_ECoG(i1(1:16))

[b2, i2] = sort(freq_region_seeg,'descend');
top_regions_seeg = unique_regions_SEEG(i1(1:15))

figure(1);clf
subplot(1,2,1)
bar(b1(2:16))
subplot(1,2,2)
bar(b2(1:15))
%% Figure 2B: derive distance - connectivity relationship in good outcome non-resected ecog vs seeg
% need to fix resected electrode indices and get ROIs from atlas to put in
% base files
for k = 1:5
ecog_dists = [];
ecog_conns = [];

seeg_dists = [];
seeg_conns = [];

for s = 1:8
    if ecog_struct(s).outcome
        
        % get patient distances
        patient_dist = ecog_distance(s).data;
        
        % delete out resected electrodes
        patient_dist(ecog_struct(s).res_elec_inds,:) = [];
        patient_dist(:,ecog_struct(s).res_elec_inds) = [];
        
        % delete out resected electrodes
        patient_conn = ecog_connectivity(s).data;
        patient_conn(ecog_struct(s).res_elec_inds,:,:) = [];
        patient_conn(:,ecog_struct(s).res_elec_inds,:) = [];
        
        % delete out self conns
        all_dist = patient_dist(:);
        self_conns = find(all_dist==0);
        all_dist(self_conns) = [];
        
        ecog_dists = [ecog_dists;all_dist];
        

            freq_conn = patient_conn(:,:,k);
            freq_conn = freq_conn(:);
            freq_conn(self_conns) = [];
            
            ecog_conns = [ecog_conns; freq_conn];

    end
    
    
    
    if seeg_struct(s).outcome
        % get patient distances
        patient_dist = seeg_distance(s).data;
        
        % delete out resected electrodes
        patient_dist(seeg_struct(s).res_elec_inds,:) = [];
        patient_dist(:,seeg_struct(s).res_elec_inds) = [];
        
        % delete out resected electrodes
        patient_conn =seeg_connectivity(s).data;
        patient_conn(seeg_struct(s).res_elec_inds,:,:) = [];
        patient_conn(:,seeg_struct(s).res_elec_inds,:) = [];
        
        % delete out self conns
        all_dist = patient_dist(:);
        self_conns = find(all_dist==0);
        all_dist(self_conns) = [];
        
        seeg_dists = [seeg_dists;all_dist];
        

            freq_conn = patient_conn(:,:,k);
            freq_conn = freq_conn(:);
            freq_conn(self_conns) = [];
            
            seeg_conns = [seeg_conns; freq_conn];

    end
    
    
end

% create figure
figure(k);clf;
subplot(1,2,1)
plot(ecog_dists,ecog_conns,'k.')
axis([0,200,0,1])
subplot(1,2,2)
plot(seeg_dists,seeg_conns,'k.')
axis([0,200,0,1])

all_dist = [ecog_dists;seeg_dists];
all_conns = [ecog_conns;seeg_conns];

[fitresult(k).data] = fit(all_dist,all_conns,'exp2');

end
%% Figure 2C: Connectivity of individual nodes

for k = 1:5
    all_ecog_conn_var = [];
    all_seeg_conn_var = [];
    
    all_ecog_nodestr = [];
    all_seeg_nodestr = [];
    
    for s = 1:8
        ecog_conn = ecog_struct(s).ii(k).data;
        seeg_conn = seeg_struct(s).ii(k).data;
        
        % variance in connectivity of individual nodes
        ecog_conn_var(s).data(k,:) = var(ecog_conn);
        seeg_conn_var(s).data(k,:) = var(seeg_conn);
        
        % calculate node strength
        ecog_nodestr(s).data(k,:) = sum(ecog_conn);
        seeg_nodestr(s).data(k,:) = sum(seeg_conn);
        
        all_ecog_conn_var = [all_ecog_conn_var, var(ecog_conn)];
        all_seeg_conn_var = [all_seeg_conn_var, var(seeg_conn)];
        
        all_ecog_nodestr = [all_ecog_nodestr, sum(ecog_conn)];
        all_seeg_nodestr = [all_seeg_nodestr, sum(seeg_conn)];
    end
    figure(k);clf;
    subplot(1,2,1)
    histogram(all_ecog_conn_var)
    xlabel('Variance in nodal connectivity')
    ylabel('Number of nodes')
    title('Group - level ECoG')
    subplot(1,2,2)
    histogram(all_seeg_conn_var)
    xlabel('Variance in nodal connectivity')
    ylabel('Number of nodes')
    title('Group - level SEEG')
end 



%% Figure 2D: Heterogeneity in node strengths

all_var_nodestr = [];

for k = 1:5
    all_nodestr_var_ecog = [];
    all_nodestr_var_seeg = [];
    
    for s = 1:8
        pt_nodestr_var_ecog = var(ecog_nodestr(s).data(k,:));
        pt_nodestr_var_seeg = var(seeg_nodestr(s).data(k,:));
        
        all_nodestr_var_ecog = [all_nodestr_var_ecog, pt_nodestr_var_ecog];
        all_nodestr_var_seeg = [all_nodestr_var_seeg, pt_nodestr_var_seeg];
    end
    
    all_var_nodestr = [all_var_nodestr, all_nodestr_var_ecog, all_nodestr_var_seeg];
    
    [p1(k),h1(k)] = ranksum(all_nodestr_var_ecog,all_nodestr_var_seeg);
    
    %figure(k);clf;
    %scatter([ones(1,8),2*ones(1,8)],[all_nodestr_var_ecog,all_nodestr_var_seeg],75,'filled')
    
end

% make figure
figure(1);clf;
x_axis = [0.5*ones(1,8), ones(1,8), 2*ones(1,8), 2.5*ones(1,8), 3.5*ones(1,8), 4*ones(1,8), 5*ones(1,8), 5.5*ones(1,8), 6.5*ones(1,8), 7*ones(1,8)];
C = zeros(80,3)+[78 172 91]/255;
C([1:8,17:24,33:40,49:56,65:72],:) = C([1:8,17:24,33:40,49:56,65:72],:)-[78 172 91]/255+[103 55 155]/255;
scatter(x_axis, all_var_nodestr,75,C,'filled')
ylabel('Variance in node strength')
%% Figure 3A: Modularity renderings
for k = 1:5
    for s = 1:8
        % ecog
        [ecog_S,ecog_Q] = modularity_und(ecog_struct(s).ii(k).data,1);
        ecog_part_coef(s).freq(k).data = participation_coef(ecog_struct(s).ii(k).data,ecog_S);
        
        ecog_modularity(s,k) = ecog_Q;
        ecog_modules(s).data(k,:) = ecog_S;
        
        % seeg
        [seeg_S,seeg_Q] = modularity_und(seeg_struct(s).ii(k).data,1);
        seeg_part_coef(s).freq(k).data = participation_coef(seeg_struct(s).ii(k).data,seeg_S);
        
        seeg_modularity(s,k) = seeg_Q;
        seeg_modules(s).data(k,:) = seeg_S;
    end
end

%% Figure 3B: Frequency dependent modularity quantification

% make figure
figure(1);clf
subplot(1,2,1)
hold on
for s = 1:8
    plot(ecog_modularity(s,:))
end
hold off

subplot(1,2,2)
hold on
for s = 1:8
    plot(seeg_modularity(s,:))
end
hold off

for k = 1:5
    mean(ecog_modularity(:,k))
    mean(seeg_modularity(:,k))
    [p2(k), h2(k)] = ranksum(ecog_modularity(:,k),seeg_modularity(:,k))
end

% group figure
figure(2);clf
all_modularity = [ecog_modularity(1:8),seeg_modularity(1:8),ecog_modularity(9:16),seeg_modularity(9:16),ecog_modularity(17:24),seeg_modularity(17:24),ecog_modularity(25:32),seeg_modularity(25:32),ecog_modularity(33:40),seeg_modularity(33:40)];
scatter(x_axis, all_modularity,75,C,'filled')
axis([0, 7.5, 0, 0.15])

%% Figure 3C: Modularity vs anatomic localization

% dice coefficient between module and AAL116 atlas that overlaps most?

% for k = 1:5
%     for s = 1:8
%         % for each module -> find to what extent it overlaps with most
%         % overlapping brain ROI
%         unique_ecog_modules = unique(ecog_modules(s).data(k,:));
%         for i = 1:length(unique_ecog_modules)
%             constituent_nodes = find(ecog_modules(s).data(k,:)==i);
%             
%             % get anatomic ROI
%             module_roi = ecog_struct(s).elec_roi(constituent_nodes);
%             unique_module_roi = length(unique(module_roi));
%             
%             num_ecog_mod_roi(s).freq(k).data(i) = unique_module_roi;
%             
%         end
%         
%         % Do with SEEG
%         unique_seeg_modules = unique(seeg_modules(s).data(k,:));
%         for i = 1:length(unique_seeg_modules)
%             constituent_nodes = find(seeg_modules(s).data(k,:)==i);
%             
%             % get anatomic ROI
%             module_roi = seeg_struct(s).elec_roi(constituent_nodes);
%             unique_module_roi = length(unique(module_roi));
%             
%             num_seeg_mod_roi(s).freq(k).data(i) = unique_module_roi;
%             
%         end
%     end
% end
% 
% % make 5 plots. X is frequency band Y is 
% 
% 
% for k = 1:5
%     all_ecog_mod_roi = [];
%     all_seeg_mod_roi = [];
%     
%     for s = 1:8
%         all_ecog_mod_roi = [all_ecog_mod_roi, num_ecog_mod_roi(s).freq(k).data];
%         all_seeg_mod_roi = [all_seeg_mod_roi, num_seeg_mod_roi(s).freq(k).data];
%         
%     end
%     
%     
%     
%     mean(all_ecog_mod_roi)
%     mean(all_seeg_mod_roi)
%     
%     figure(k);clf;
%     subplot(1,2,1)
%     histogram(all_ecog_mod_roi,8)
%     subplot(1,2,2)
%     histogram(all_seeg_mod_roi,8)
% end


%figure(6);clf
%scatter(x_axis, all_mod_roi,75,C,'filled')


%% Do localization based on Yeo atlas
% Hypothesize that modules fall into single system more commonly on SEEG
% than ECoG. This would imply significant amount of connectivity linkage is
% in areas known to be linked that may not be relevant to epilepsy. 

% region
% fileID = fopen('/Users/jbernabei/Downloads/Atlas_labels/Yeo2011_7Networks_ColorLUT.txt');
% atlas_info1 = textscan(fileID,'%s %s %d');
% all_inds1 = [double(atlas_info1{3})];
% all_locs1 = [atlas_info1{2}];
% 
% for s = 1:8
%     clear atlas_yeo_ecog
%     clear atlas_yeo_seeg
%     
%    % Do electrode localization for atlases based on MNI coordinates
%    [mni_coords1, ecog_systems, NN_flag1] = nifti_values(ecog_struct(s).mni_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'); 
%    % Do electrode localization for atlases based on MNI coordinates
%    [mni_coords2, seeg_systems, NN_flag2] = nifti_values(seeg_struct(s).mni_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'); 
% 
%     % ecog
%     num_elecs = size(mni_coords1,1);
%     for i  = 1:num_elecs
%         yeo_string = find(ecog_systems(i)==all_inds1);
%         if isempty(yeo_string) 
%             atlas_yeo_ecog{i,1} = 'n/a';
%         else
%             atlas_yeo_ecog{i,1} = all_locs1{yeo_string};
%         end
%     end
%     
%     % seeg
%     num_elecs = size(mni_coords2,1);
%     for i  = 1:num_elecs
%         yeo_string = find(seeg_systems(i)==all_inds1);
%         if isempty(yeo_string) 
%             atlas_yeo_seeg{i,1} = 'n/a';
%         else
%             atlas_yeo_seeg{i,1} = all_locs1{yeo_string};
%         end
%     end
%     
%     % now we want to find dice of each module with what system it primarily
%     % represents
%     for k = 1:5
%         
%         % find unique ecog modules
%         unique_ecog_modules = unique(ecog_modules(s).data(k,:));
%         for i = 1:length(unique_ecog_modules)
%             constituent_nodes = find(ecog_modules(s).data(k,:)==i);
%             
%             for r = 1:7
%                 system_nodes = find(r==ecog_systems);
%                 ecog_system_overlap(s).freq(k).data(i,r) = dice(constituent_nodes,system_nodes);
%             end
%             
%         end
%         
%         % Do with SEEG
%         unique_seeg_modules = unique(seeg_modules(s).data(k,:));
%         for i = 1:length(unique_seeg_modules)
%             constituent_nodes = find(seeg_modules(s).data(k,:)==i);
%             
%             for r = 1:7
%                 system_nodes = find(r==seeg_systems);
%                 seeg_system_overlap(s).freq(k).data(i,r) = dice(constituent_nodes,system_nodes);
%             end
%             
%         end
%         
%         max_overlap_ecog(s,k) = mean(max(ecog_system_overlap(s).freq(k).data'));
%         max_overlap_seeg(s,k) = mean(max(seeg_system_overlap(s).freq(k).data'));
%         
%     end
%     
% end


%% Figure 4B: Node strength localization

% figure out avg number of ecog elecs resected
for s = 1:8
    num_resect_ecog(s) = length(ecog_struct(s).res_elec_inds);
    num_resect_seeg(s) = length(seeg_struct(s).res_elec_inds);
end

mean(num_resect_ecog)
mean(num_resect_seeg)

all_all_loc_qual = [];


for k = 1:5
    
    all_seeg_target = [];
    
    all_loc_qual_ecog = [];
    all_loc_qual_seeg = [];
    
    for s = 1:8
        W_ecog = ecog_struct(s).ii(k).data;
        W_seeg = seeg_struct(s).ii(k).data;
        
        mni_coords_ecog = ecog_struct(s).mni_coords;
        mni_coords_seeg = seeg_struct(s).mni_coords;
        
        num_targets_ecog = 14;
        num_targets_seeg = 9;
        
        [targets_ecog1, targets_ecog2] = localize_EZ_II(W_ecog, num_targets_ecog, mni_coords_ecog);
        [targets_seeg1, targets_seeg2] = localize_EZ_II(W_seeg, num_targets_seeg, mni_coords_seeg);
        
        loc_qual_ecog(s).freq(k).data = dice(targets_ecog2,ecog_struct(s).res_elec_inds);
        loc_qual_seeg(s).freq(k).data = dice(targets_seeg2,seeg_struct(s).res_elec_inds);
        
        all_loc_qual_ecog = [all_loc_qual_ecog, loc_qual_ecog(s).freq(k).data];
        all_loc_qual_seeg = [all_loc_qual_seeg, loc_qual_seeg(s).freq(k).data];
        
        seeg_target_roi(s).freq(k).data = seeg_struct(s).elec_roi(targets_seeg1);
        
        all_seeg_target = [all_seeg_target; seeg_target_roi(s).freq(k).data];
        
        record_seeg_targets(s).freq(k).data = targets_seeg1;
        
        
    end
    
    all_all_loc_qual = [all_all_loc_qual, all_loc_qual_ecog, all_loc_qual_seeg];
    
%     figure(k);clf;
%     scatter([ones(1,8),2*ones(1,8)],[all_loc_qual_ecog,all_loc_qual_seeg])
    
        unique_seeg_targets = unique(all_seeg_target)
        for r = 1:length(unique_seeg_targets)
            num_target(r) = sum(strcmp(unique_seeg_targets(r),all_seeg_target));
        end
        
        num_target
        if k==4
        else
        clear num_target
        end

end


figure(1);clf
scatter(x_axis, all_all_loc_qual,75,C,'filled','jitter','on', 'jitterAmount', 0.1)
axis([0, 7.5, 0, 1.25])

%% example rendering for missed seeg localization
% blue is true resection zone
% red is incorrect proposal
% Now find top 25 for ECoG
[b1, i1] = sort(num_target,'descend');
top_targets_seeg = unique_seeg_targets(i1(1:10));
figure(1);clf
bar(b1(1:10))

    
    % Load in patient ID
    ID_seeg = seeg_struct(4).ID;
    
    % Set output path
    output_seeg = sprintf('%s/',strcat(base_path, ID_seeg));
    
    % Set filename
    seeg_filename = sprintf('localization_rendering_%s.jpg',ID_seeg);
    
    % Get mni coordinates
    mni_seeg = seeg_struct(4).mni_coords;
    
    % Get color data
    color_seeg = zeros(size(mni_seeg));
    color_seeg(seeg_struct(4).res_elec_inds) = -1;
    color_seeg(record_seeg_targets(4).freq(4).data) = 1;
    
    % Get color data
    size_seeg = ones(size(mni_seeg));
    
    % Create rendering
    create_rendering(output_seeg, seeg_filename, mni_seeg, color_seeg, size_seeg)


%% Figure 5A: Reduce networks by anatomic localization
% loop through patients
for s = 1:8
    
    % find unique brain regions
    unique_roi_ecog(s).data = unique(ecog_struct(s).elec_roi);
    unique_roi_seeg(s).data = unique(seeg_struct(s).elec_roi);
    
    a = 0;
    
    % loop through unique brain regions
    for i = 1:length(unique_roi_ecog(s).data)
        
        % find which nodes are in these
        nodes_1 = find(strcmp(ecog_struct(s).elec_roi, unique_roi_ecog(s).data(i)));

        % find centroid of these regions
        centroid_1 = mean(ecog_struct(s).mni_coords(nodes_1,:),1);
        pt_centr_ecog(s).data(i,:) = centroid_1;
        
        % check if its a resected region
        if sum(sum(nodes_1'==ecog_struct(s).res_elec_inds))>0
            a = a+1;
            resect_region_ecog(s).data(a) = i;
        end
        
        % loop through regions again
        for j = 1:length(unique_roi_ecog(s).data)
            
            % find second set of nodes
            nodes_2 = find(strcmp(ecog_struct(s).elec_roi, unique_roi_ecog(s).data(j)));
            
            % find second set of centroidss
            centroid_2 = mean(ecog_struct(s).mni_coords(nodes_2,:));
            if i==j
            else
                for k = 1:5
                    reduced_adj_ecog(s).freq(k).data(i,j) = mean(mean(ecog_struct(s).ii(k).data(nodes_1,nodes_2)));
                end 
            end
        end
    end
    
    a = 0;
    
    % loop through unique brain regions
    for i = 1:length(unique_roi_seeg(s).data)
        
        % find which nodes are in these
        nodes_1 = find(strcmp(seeg_struct(s).elec_roi, unique_roi_seeg(s).data(i)));

        % find centroid of these regions
        centroid_1 = mean(seeg_struct(s).mni_coords(nodes_1,:),1);
        pt_centr_seeg(s).data(i,:) = centroid_1;
        
        % check if its a resected region
        if sum(sum(nodes_1==seeg_struct(s).res_elec_inds))>0
            a = a+1;
            resect_region_seeg(s).data(a) = i;
        end
        
        % loop through regions again
        for j = 1:length(unique_roi_seeg(s).data)
            
            % find second set of nodes
            nodes_2 = find(strcmp(seeg_struct(s).elec_roi, unique_roi_seeg(s).data(j)));
            
            % find second set of centroidss
            centroid_2 = mean(seeg_struct(s).mni_coords(nodes_2,:));
            if i==j
            else
                for k = 1:5
                    reduced_adj_seeg(s).freq(k).data(i,j) = mean(mean(seeg_struct(s).ii(k).data(nodes_1,nodes_2)));
                end 
            end
        end
    end
    
    
end

%% Figure 5 reduced renderings

% Loop through all patients
for s = 1:8
    
    % Load in patient ID
    ID_ecog = ecog_struct(s).ID;
    ID_seeg = seeg_struct(s).ID;
    
    % Set output path
    output_ecog = sprintf('%s/',strcat(base_path, ID_ecog));
    output_seeg = sprintf('%s/',strcat(base_path, ID_seeg));
    
    % Set filename
    ecog_filename = sprintf('reduced_rendering_%s.jpg',ID_ecog);
    seeg_filename = sprintf('reduced_rendering_%s.jpg',ID_seeg);
    
    % Get mni coordinates
    mni_ecog = pt_centr_ecog(s).data;
    mni_seeg = pt_centr_seeg(s).data;
    
    % Get color data
    color_ecog = -1*ones(size(mni_ecog));
    color_seeg = -1*ones(size(mni_seeg));
    
    % Get color data
    size_ecog = ones(size(mni_ecog));
    size_seeg = ones(size(mni_seeg));
    
    % Create rendering
    create_rendering(output_ecog, ecog_filename, mni_ecog, color_ecog, size_ecog)
    create_rendering(output_seeg, seeg_filename, mni_seeg, color_seeg, size_seeg)
    
end

%% Figure 5B: Compare distancce

all_ecog_reduced_dist = [];
all_seeg_reduced_dist = [];

for s = 1:8
    num_roi_ecog = size(pt_centr_ecog(s).data,1);
    num_roi_seeg = size(pt_centr_seeg(s).data,1);

    % loop through ecog
    for i = 1:num_roi_ecog
        for j = 1:num_roi_ecog
            roi1 = pt_centr_ecog(s).data(i,:);
            roi2 = pt_centr_ecog(s).data(j,:);
            
            dist_ecog_reduced(s).data(i,j) = sqrt(sum((roi1-roi2).^2));
        end
    end

    % loop through seeg
    for i = 1:num_roi_seeg
        for j = 1:num_roi_seeg
            roi1 = pt_centr_seeg(s).data(i,:);
            roi2 = pt_centr_seeg(s).data(j,:);
            
            dist_seeg_reduced(s).data(i,j) = sqrt(sum((roi1-roi2).^2));
        end
    end
    
    all_ecog_reduced_dist = [all_ecog_reduced_dist; dist_ecog_reduced(s).data(:)];
    all_seeg_reduced_dist = [all_seeg_reduced_dist; dist_seeg_reduced(s).data(:)];

end

all_ecog_reduced_dist(find(all_ecog_reduced_dist==0)) = [];
all_seeg_reduced_dist(find(all_seeg_reduced_dist==0)) = [];

mean(all_ecog_reduced_dist)
skewness(all_ecog_reduced_dist)

mean(all_seeg_reduced_dist)
skewness(all_seeg_reduced_dist)

% make figure
figure(1);clf
subplot(1,2,1)
histogram(all_ecog_reduced_dist)
xlabel('Internodal distance (mm)')
ylabel('Number of edges')
title('Group - level reduced ECoG')
subplot(1,2,2)
histogram(all_seeg_reduced_dist)
xlabel('Internodal distance (mm)')
ylabel('Number of edges')
title('Group - level reduced SEEG')

%% Figure 5C: Compare connectivity heterogeneity

for k = 1:5
    connvar_reduced_ecog = [];
    connvar_reduced_seeg = [];
    
    for s = 1:8
        connvar_reduced_ecog = [connvar_reduced_ecog, var(reduced_adj_ecog(s).freq(k).data)];
        connvar_reduced_seeg = [connvar_reduced_seeg, var(reduced_adj_seeg(s).freq(k).data)];
    end
    
    % make figure
    figure(k);clf
    subplot(1,2,1)
    histogram(connvar_reduced_ecog)
    xlabel('Variance in nodal connectivity')
    ylabel('Number of nodes')
    title('Group - level reduced ECoG')
    subplot(1,2,2)
    histogram(connvar_reduced_seeg)
    xlabel('Variance in nodal connectivity')
    ylabel('Number of nodes')
    title('Group - level reduced SEEG')
    
end

%% Figure 5D: Compare nodal heterogeneity

all_red_nodevar = [];

for k = 1:5
    all_nodestr_var_ecog = [];
    all_nodestr_var_seeg = [];
    
    for s = 1:8
        
        ecog_reduced_nodestr(s).data(k,:) = sum(reduced_adj_ecog(s).freq(k).data);
        seeg_reduced_nodestr(s).data(k,:) = sum(reduced_adj_seeg(s).freq(k).data);
        
        pt_nodestr_var_ecog = var(ecog_reduced_nodestr(s).data(k,:))
        pt_nodestr_var_seeg = var(seeg_reduced_nodestr(s).data(k,:))
        
        all_nodestr_var_ecog = [all_nodestr_var_ecog, pt_nodestr_var_ecog];
        all_nodestr_var_seeg = [all_nodestr_var_seeg, pt_nodestr_var_seeg];
    end
    
    all_red_nodevar = [all_red_nodevar, all_nodestr_var_ecog, all_nodestr_var_seeg];
    
    %figure(k);clf;
    %scatter([ones(1,8),2*ones(1,8)],[all_nodestr_var_ecog,all_nodestr_var_seeg])
end

figure(1);clf
scatter(x_axis, all_red_nodevar,75,C,'filled','jitter','on', 'jitterAmount', 0.1)
axis([0, 7.5, 0, 6])

%% Figure 5E: Compare modularity

for k = 1:5
    for s = 1:8
        [ecog_rS,ecog_rQ] = modularity_und(reduced_adj_ecog(s).freq(k).data);
        ecog_r_part_coef(s).freq(k).data = participation_coef(reduced_adj_ecog(s).freq(k).data,ecog_rS);

        ecog_r_modularity(s,k) = ecog_rQ;
        ecog_r_modules(s).data(k,:) = ecog_rS;

        % seeg
        [seeg_rS,seeg_rQ] = modularity_und(reduced_adj_seeg(s).freq(k).data,1);
        seeg_r_part_coef(s).freq(k).data = participation_coef(reduced_adj_seeg(s).freq(k).data,seeg_rS);

        seeg_r_modularity(s,k) = seeg_rQ;
        seeg_r_modules(s).data(k,:) = seeg_rS;
    end
end

for k = 1:5
    mean(ecog_r_modularity(:,k))
    mean(seeg_r_modularity(:,k))
    [p2(k), h2(k)] = ranksum(ecog_r_modularity(:,k),seeg_r_modularity(:,k))
end

% make figure
figure(1);clf
all_modularity = [ecog_r_modularity(1:8),seeg_r_modularity(1:8),ecog_r_modularity(9:16),seeg_r_modularity(9:16),ecog_r_modularity(17:24),seeg_r_modularity(17:24),ecog_r_modularity(25:32),seeg_r_modularity(25:32),ecog_r_modularity(33:40),seeg_r_modularity(33:40)];
scatter(x_axis, all_modularity,75,C,'filled')
axis([0, 7.5, 0, 0.15])

%% Figure 6: Compare localization quality

all_all_loc_qual_r = [];

for k = 1:5
    
    all_seeg_r_target = [];
    
    all_loc_qual_ecog = [];
    all_loc_qual_seeg = [];
    
    adjust = @(x) fitresult(k).data.a.*exp(x.*fitresult(k).data.b)+fitresult(k).data.c.*exp(x.*fitresult(k).data.d);
    
    for s = 1:8
        W_ecog = reduced_adj_ecog(s).freq(k).data;
        W_seeg = reduced_adj_seeg(s).freq(k).data;
        
        dist_matrix_ecog = dist_ecog_reduced(s).data;
        dist_matrix_seeg = dist_seeg_reduced(s).data;
        
        normal_ecog = adjust(dist_matrix_ecog);
        normal_seeg = adjust(dist_matrix_seeg);
        
        W_ecog = W_ecog;%-normal_ecog;
        W_seeg = W_seeg;%-normal_seeg;
        
        delete_ind_ecog = [];
        delete_ind_seeg = [];
        for q = 1:length(unique_roi_ecog(s).data)
            if strcmp(unique_roi_ecog(s).data{q},'White_Matter')||strcmp(unique_roi_ecog(s).data{q},'n/a')
                delete_ind_ecog = [delete_ind_ecog, q];
            end
        end
        
        for q = 1:length(unique_roi_seeg(s).data)
            if strcmp(unique_roi_seeg(s).data{q},'White_Matter')||strcmp(unique_roi_seeg(s).data{q},'n/a')
                delete_ind_seeg = [delete_ind_seeg, q];
            end
        end
        
        mni_coords_ecog = pt_centr_ecog(s).data;
        mni_coords_seeg = pt_centr_seeg(s).data;
        
        %W_ecog(delete_ind_ecog,:) = [];
        %W_ecog(:,delete_ind_ecog) = [];
        %mni_coords_ecog(delete_ind_ecog,:) = [];
        
        W_seeg(delete_ind_seeg,:) = [];
        W_seeg(:,delete_ind_seeg) = [];
        mni_coords_seeg(delete_ind_seeg,:) = [];
        
        num_targets_ecog = 5;
        num_targets_seeg = 3;
        
        [targets_r_ecog1, targets_r_ecog2] = localize_EZ_II(W_ecog, num_targets_ecog, mni_coords_ecog);
        [targets_r_seeg1, targets_r_seeg2] = localize_EZ_II(W_seeg, num_targets_seeg, mni_coords_seeg);
        
        loc_qual_ecog(s).freq(k).data = dice(targets_r_ecog2,resect_region_ecog(s).data);
        loc_qual_seeg(s).freq(k).data = dice(targets_r_seeg2,resect_region_seeg(s).data);
        
        all_loc_qual_ecog = [all_loc_qual_ecog, loc_qual_ecog(s).freq(k).data];
        all_loc_qual_seeg = [all_loc_qual_seeg, loc_qual_seeg(s).freq(k).data];
        
        all_seeg_r_target = [all_seeg_r_target; unique_roi_seeg(s).data(targets_r_seeg1)];
    end
    
    all_loc_qual_ecog;
    all_loc_qual_seeg;
    
    all_all_loc_qual_r = [all_all_loc_qual_r, all_loc_qual_ecog, all_loc_qual_seeg];
    
    
        unique_seeg_r_targets = unique(all_seeg_r_target)
        for r = 1:length(unique_seeg_r_targets)
            num_target_r(r) = sum(strcmp(unique_seeg_r_targets(r),all_seeg_r_target));
        end
        
        num_target_r
        clear num_target_r
    
    %figure(k);clf;
    %scatter([ones(1,8),2*ones(1,8)],[all_loc_qual_ecog,all_loc_qual_seeg])
    
end

% make figure
figure(1);clf
scatter(x_axis, all_all_loc_qual_r,75,C,'filled','d','jitter','on', 'jitterAmount', 0.1)
axis([0, 7.5, 0, 1.25])

%%

%% Correct adjacency matrices for distance. 

for k = 1:5
    
    adjust = @(x) fitresult(k).data.a.*exp(x.*fitresult(k).data.b)+fitresult(k).data.c.*exp(x.*fitresult(k).data.d);
    
    all_loc_c_qual_ecog = [];
    all_loc_c_qual_seeg = [];
    for s = 1:8
        W_ecog = ecog_struct(s).ii(k).data;
        W_seeg = seeg_struct(s).ii(k).data;
        
        mni_coords_ecog = ecog_struct(s).mni_coords;
        mni_coords_seeg = seeg_struct(s).mni_coords;
        
        num_targets_ecog = 5;
        num_targets_seeg = 3;
        
        normal_ecog = alphatheta_adjust(ecog_distance(s).data);
        normal_seeg = alphatheta_adjust(seeg_distance(s).data);
        
        corrected_ecog = W_ecog;%-normal_ecog;
        corrected_seeg = W_seeg;%-normal_seeg;
        
        [targets_ecog1, targets_ecog2] = localize_EZ_II(corrected_ecog, num_targets_ecog, mni_coords_ecog);
        [targets_seeg1, targets_seeg2] = localize_EZ_II(corrected_seeg, num_targets_seeg, mni_coords_seeg);
        
        loc_qual_c_ecog(s).freq(k).data = dice(targets_ecog1,ecog_struct(s).res_elec_inds);
        loc_qual_c_seeg(s).freq(k).data = dice(targets_seeg1,seeg_struct(s).res_elec_inds);
        
        all_loc_c_qual_ecog = [all_loc_c_qual_ecog, loc_qual_c_ecog(s).freq(k).data];
        all_loc_c_qual_seeg = [all_loc_c_qual_seeg, loc_qual_c_seeg(s).freq(k).data];
    end
    all_loc_c_qual_ecog
    all_loc_c_qual_seeg
end
