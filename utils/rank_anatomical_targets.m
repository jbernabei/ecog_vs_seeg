function [top_regions, plot_data, total_elecs, all_loc_matrix, bilat_score,ipsi_score] = rank_anatomical_targets(patientID, patient_roi, laterality, targets, atlas_inds, atlas_locs, lobe_table)
    
    % get number of patients
    num_patients = length(patientID);
    
    % loop through each patient
    for pt = 1:num_patients
        clear loc_table
        clear loc_matrix
        
        % 
        [I,J] = find(atlas_inds==patient_roi{pt});
        patient_regions = {atlas_locs{I}}';
       
        % correct R/L into ipsilateral / contralateral
        if strcmp(laterality{pt},'R')
            ipsilateral = contains(patient_regions,'_R');
            contralateral = contains(patient_regions,'_L');
        else
            ipsilateral = contains(patient_regions,'_L');
            contralateral = contains(patient_regions,'_R');
        end
        
        if strcmp(targets{pt},'Temporal') && strcmp(laterality{pt},'R')
            res_region = 1;
        elseif strcmp(targets{pt},'Temporal') && strcmp(laterality{pt},'L')
            res_region = 2;
        elseif strcmp(targets{pt},'Frontal') && strcmp(laterality{pt},'R')
            res_region = 3;
        elseif strcmp(targets{pt},'Frontal') && strcmp(laterality{pt},'L')
            res_region = 4;
        elseif strcmp(targets{pt},'Parietal') && strcmp(laterality{pt},'R')
            res_region = 5;
        elseif strcmp(targets{pt},'Parietal') && strcmp(laterality{pt},'L')
            res_region = 6;
        elseif strcmp(targets{pt},'Occipital') && strcmp(laterality{pt},'R')
            res_region = 7;
        elseif strcmp(targets{pt},'Occipital') && strcmp(laterality{pt},'L')
            res_region = 8;
        elseif strcmp(targets{pt},'Insular') && strcmp(laterality{pt},'R')
            res_region = 9;
        elseif strcmp(targets{pt},'Insular') && strcmp(laterality{pt},'L')
            res_region = 10;
        else
            res_region = 0;
        end
     
        %
        for r = 1:length(patient_regions)
            processed_roi = 'n/a';

            roi_name = patient_regions{r};
            
            % we want to make plots of proportion of interhemispheric edges,
            % interlobar edges, and intralobar edges
        
            % step 1: use atlas table to map to R/L_lobe
            % use key codes - 1/2/3/4/5/6/7/8
            loc_node = lobe_table{find(strcmp(lobe_table{:,1},roi_name)),3};
        
            if isempty(loc_node)
                loc_node = 'Empty';
            end
            
            if strcmp(loc_node,'R_Temporal') || strcmp(loc_node,'R_MTL')
                loc_table(r,1) = 1;
            elseif strcmp(loc_node,'L_Temporal') || strcmp(loc_node,'L_MTL')
                loc_table(r,1) = 2;
            elseif strcmp(loc_node,'R_Frontal') || strcmp(loc_node,'R_MFL')
                loc_table(r,1) = 3;
            elseif strcmp(loc_node,'L_Frontal') || strcmp(loc_node,'L_MFL')
                loc_table(r,1) = 4;
            elseif strcmp(loc_node,'R_Parietal')
                loc_table(r,1) = 5;
            elseif strcmp(loc_node,'L_Parietal')
                loc_table(r,1) = 6;
            elseif strcmp(loc_node,'R_Occipital')
                loc_table(r,1) = 7;
            elseif strcmp(loc_node,'L_Occipital')
                loc_table(r,1) = 8;
            elseif strcmp(loc_node,'R_Insular')
                loc_table(r,1) = 9;
            elseif strcmp(loc_node,'L_Insular')
                loc_table(r,1) = 10;
            elseif strcmp(roi_name,'White_Matter')
                loc_table(r,1) = 0.5;
            end

            new_roi = strip(roi_name,'right','R');
            new_roi = strip(new_roi,'right','L');

            if ipsilateral(r)==1
                processed_roi = strcat(new_roi,'ipsilateral');
            elseif contralateral(r)==1
                processed_roi = strcat(new_roi,'contralateral');
            end

            all_processed_roi(pt).data{r} = processed_roi;

        end
       
        % create a bilaterality score: fraction of nodes contra/ipsilateral
        num_ipsi = sum(contains(all_processed_roi(pt).data,'ipsilateral'));
        num_cont = sum(contains(all_processed_roi(pt).data,'contralateral'));
        bilat_score(pt) = min([num_ipsi,num_cont])./max([num_ipsi,num_cont]);
        
        % check for intralobar, interlobar, interhemispheric
        for r1 = 1:length(patient_regions)
            for r2 = 1:length(patient_regions)
                this_loc_1 = loc_table(r1);
                this_loc_2 = loc_table(r2);
                
                if this_loc_1==this_loc_2
                    loc_matrix(r1,r2) = 1; % intralobar
                elseif mod(this_loc_1,2)==mod(this_loc_2,2)
                    loc_matrix(r1,r2) = 2; % interlobar
                elseif (mod(this_loc_1,2) == 1) && (mod(this_loc_2,2) == 0) || (mod(this_loc_1,2) == 0) && (mod(this_loc_2,2) == 1)
                    loc_matrix(r1,r2) = 3; % interhemispheric
                else
                    loc_matrix(r1,r2) = 4; % white matter
                end
            end 
        end
        
    all_loc_matrix(pt).data = loc_matrix;
    
    % create a ipsilateral focality score: fraction of ipsilateral nodes in
    % lobe with greatest number of resected electrodes
    ipsi_inds = contains(all_processed_roi(pt).data,'ipsilateral');
    ipsi_score(pt) = sum(loc_table(ipsi_inds)==res_region)./sum(ipsi_inds);
    
    % <or we can just define lobe of primary target>
    
    
    end
    
    
    
    % Now we have to go through all patients and determine ranked list of ROI
    % overall in ECoG and apply that list to SEEG
    all_roi = [];

    for pt = 1:num_patients
        all_roi = [all_roi; all_processed_roi(pt).data'];
    end
    
    num_n_a = find(strcmp(all_roi,'n/a'));
    all_roi(num_n_a) = [];

    unique_regions = unique(all_roi)

    for r = 1:length(unique_regions)
        freq_region(r) = sum(strcmp(all_roi,unique_regions{r}));
    end

    % Now find top 15 for ECoG
    [b1, i1] = sort(freq_region,'descend');
    top_regions = unique_regions(i1(1:15));
    plot_data = b1(1:15);
    
    total_elecs = sum(b1);

end
    
    