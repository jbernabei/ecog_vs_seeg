% ecog vs seeg analysis

% create structure for each patient with the following information:
%   1. pt name
%   2. pt outcome (most recent)
%   3. adj matrices for one seizure of each type in each freq. band
%   4. patient electrode names
%   5. patient resected electrode names
%   6. patient resected electrode indices
%   7. patient electrode mni coordinates
%   8. patient electrode ROI
%   9. interictal adj matrix from 1 hour activity average

patient_name = 'HUP078';
patient_outcome = 0;

path_to_adj = '/Users/jbernabei/Documents/PhD_Research/ecog_vs_seeg/HUP078'
path_to_coords = '/Users/jbernabei/Documents/PhD_Research/ecog_vs_seeg/HUP078'
path_to_elec_labels = '/Users/jbernabei/Documents/PhD_Research/Virtual_Resection'

ignore_electrodes = {'EKG1','EKG2','AMY1','AMY2','HIPP1','HIPP2','HIPP3','HIPP4','LG1'}
resected_electrodes = {'/Users/jbernabei/Documents/PhD_Research/Virtual_Resection/ECoG_res_elecs'}

patient_struct.ID = patient_name
patient_struct.outcome = patient_outcome

%% get adj matrices
file_names_list = dir(path_to_adj);
num_files = length(file_names_list);
a = 0;
for f = 1:num_files
    file_name = file_names_list(f).name;
    if strfind(file_name,'mat.multiband.npz')
        a = a+1;
        unzip(sprintf('%s/%s',path_to_adj,file_name));
        patient_adj(1).freq(1).data = readNPY('all_adj_beta.npy');
        patient_adj(1).freq(2).data = readNPY('all_adj_lowgamma.npy');
        patient_adj(1).freq(3).data = readNPY('all_adj_highgamma.npy');
        patient_adj(1).freq(4).data = readNPY('all_adj_broadband_CC.npy');

        delete all_adj_beta.npy all_adj_lowgamma.npy all_adj_highgamma.npy all_adj_broadband_CC.npy 

    elseif strfind(file_name,'1.multiband.npz')
        unzip(sprintf('%s/%s',path_to_adj,file_name));
        patient_adj(2).freq(1).data = readNPY('all_adj_beta.npy');
        patient_adj(2).freq(2).data = readNPY('all_adj_lowgamma.npy');
        patient_adj(2).freq(3).data = readNPY('all_adj_highgamma.npy');
        patient_adj(2).freq(4).data = readNPY('all_adj_broadband_CC.npy');
    end
end

%%

for k = 1:4

    patient_struct.II.mean.freq(k).data = nanmean(patient_adj(1).freq(k).data,3);
    patient_struct.II.var.freq(k).data = nanvar(patient_adj(1).freq(k).data,[],3);
    
    patient_struct.seizure.freq(k).data = patient_adj(1).freq(k).data;

end

%% 106 channels

elec_labels_raw = readtable('/Users/jbernabei/Documents/PhD_Research/Virtual_Resection/HUP078_electrode_labels.csv')
elec_labels_raw = elec_labels_raw{:,5}

ignore_inds = [];

for i = 1:length(elec_labels_raw)
    if sum(strcmp(elec_labels_raw{i},ignore_electrodes))>0
        ignore_inds = [ignore_inds, i]
    end
end

elec_labels = elec_labels_raw;
elec_labels(ignore_inds) = [];

patient_struct.elec_labels = elec_labels;

%% resection zone

resected_elecs = readtable('/Users/jbernabei/Documents/PhD_Research/Virtual_Resection/ECoG_res_elecs/HUP078_resected_electrodes_0.csv');

patient_struct.resected_elecs = resected_elecs{:,2};
patient_struct.res_elec_inds = resected_elecs{:,1}+1;

%% coordinates
mni_coords_raw = readtable(sprintf('%s/electrodenames_coordinates_mni.csv',path_to_coords));
t1_coords_raw = readtable(sprintf('%s/electrodenames_coordinates_native_and_T1.csv',path_to_coords));

for e = 1:size(t1_coords_raw,1)
    elec_name = t1_coords_raw{e,1};
    elec_ind = find(strcmp(patient_struct.elec_labels,elec_name));
    elec_roi(elec_ind,1) = t1_coords_raw{e,2};
    if ~isempty(elec_ind)
    elec_mni_coords(elec_ind,:) = mni_coords_raw{e,2:4}; 
    end
end

patient_struct.mni_coords = elec_mni_coords;
patient_struct.elec_roi = elec_roi;

%% make rendering
 new_coords = patient_struct.mni_coords;
    
    rendering_filename = sprintf('%s_mni.jpg',patient_name);

    for ch = 1:size(new_coords,1)

        % Get coordinate data
        subject_data(ch,1:3) = new_coords(ch,:);

        subject_data(ch,4) = -1; % bad channel is red

        subject_data(ch,5) = 1; 
    end

    % write the node file
    dlmwrite(sprintf('mni_reg_%s.node',patient_name),subject_data,'delimiter',' ','precision',5)

    % make the render
        BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',sprintf('mni_reg_%s.node',patient_name),'check_reg.mat',rendering_filename)

%%
simulated_interictal(1).data = mean(patient_struct(18).freq(1).data(:,:,1:end/4),3);
simulated_interictal(2).data = mean(patient_struct(18).freq(2).data(:,:,1:end/4),3);
simulated_interictal(3).data = mean(patient_struct(18).freq(3).data(:,:,1:end/4),3);
simulated_interictal(4).data = mean(patient_struct(18).freq(4).data(:,:,1:end/4),3);
%%
figure(2);clf
subplot(2,2,1)
imagesc(simulated_interictal(1).data)
subplot(2,2,2)
imagesc(simulated_interictal(2).data)
subplot(2,2,3)
imagesc(simulated_interictal(3).data)
subplot(2,2,4)
imagesc(simulated_interictal(4).data)
%%
figure(2);clf
subplot(2,2,1)
histogram(simulated_interictal(1).data(:))
axis([0.05 0.8, 0, 3000])
subplot(2,2,2)
histogram(simulated_interictal(2).data(:))
axis([0.05 0.8, 0, 3000])
subplot(2,2,3)
histogram(simulated_interictal(3).data(:))
axis([0.05 0.8, 0, 3000])
subplot(2,2,4)
histogram(simulated_interictal(4).data(:))
axis([0.05 0.8, 0, 3000])

%% 
a = 0;
for i = 1:106
    for j = 1:106
        a = a+1;
        elec1 = patient_struct.mni_coords(i,:);
        elec2 = patient_struct.mni_coords(j,:);
        inter_elec_distance(a) = sum(sqrt((elec1-elec2).^2));
    end
end

inter_elec_distance(inter_elec_distance==0) = []
%% 
patient_name = 'HUP133'

simulated_interictal(1).data = mean(patient_struct(6).freq(1).data(:,:,1:end/4),3);
simulated_interictal(2).data = mean(patient_struct(6).freq(2).data(:,:,1:end/4),3);
simulated_interictal(3).data = mean(patient_struct(6).freq(3).data(:,:,1:end/4),3);
simulated_interictal(4).data = mean(patient_struct(6).freq(4).data(:,:,1:end/4),3);

for i = 1:4
    [S,Q] = modularity_und(simulated_interictal(i).data,1)
    
    new_coords = patient_struct(6).elec_coords_mni(1:77,:);
    
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

