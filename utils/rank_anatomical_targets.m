function [top_regions, plot_data, total_elecs] = rank_anatomical_targets(patientID, patient_roi, laterality, atlas_inds, atlas_locs, implant_type, target_type, plot_name)
    
    % get number of patients
    num_patients = length(patientID);
    
    % loop through each patient
    for pt = 1:num_patients
        
        % 
        [I,J] = find(atlas_inds==patient_roi{pt});
        patient_regions = {atlas_locs{I}}';
    
        if strcmp(laterality{pt},'R')
            ipsilateral = contains(patient_regions,'_R');
            contralateral = contains(patient_regions,'_L');
        else
            ipsilateral = contains(patient_regions,'_L');
            contralateral = contains(patient_regions,'_R');
        end
     
        %
        for r = 1:length(patient_regions)
            processed_roi = 'n/a';

            roi_name = patient_regions{r};

            new_roi = strip(roi_name,'right','R');
            new_roi = strip(new_roi,'right','L');

            if ipsilateral(r)==1
                processed_roi = strcat(new_roi,'ipsilateral');
            elseif contralateral(r)==1
                processed_roi = strcat(new_roi,'contralateral');
            end

            all_processed_roi(pt).data{r} = processed_roi;

        end
        
    end




    % Now we have to go through all patients and determine ranked list of ROI
    % overall in ECoG and apply that list to SEEG
    all_roi = [];

    for pt = 1:num_patients
        all_roi = [all_roi; all_processed_roi(pt).data'];
    end
    
    num_n_a = find(strcmp(all_roi,'n/a'));
    all_roi(num_n_a) = [];

    unique_regions = unique(all_roi);

    for r = 1:length(unique_regions)
        freq_region(r) = sum(strcmp(all_roi,unique_regions{r}));
    end

    % Now find top 15 for ECoG
    [b1, i1] = sort(freq_region,'descend');
    top_regions = unique_regions(i1(1:15));
    plot_data = b1(1:15);
    
    total_elecs = sum(b1);
    
%     histogram_plot = bar(b1(1:15));
%     title(sprintf('Anatomical distribution of %s %s patients',implant_type,target_type));
%     set(gca,'xtick',(1:15),'xticklabel',top_regions);
%     xtickangle(45);
%     
%     saveas(gcf,sprintf('output/supplemental_figures/%s',plot_name));
%     close(gcf);

end
    
    