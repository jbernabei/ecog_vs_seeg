function [similarity_index] = network_similarity(network_modules, mni_coordinates, comparison_type)
    
    % get number of patients
    num_patients = length(mni_coordinates);
    
    for pt = 1:num_patients
        
        % extract coordinates
        mni_coords = mni_coordinates{pt};
        
        % compute systems
        [mni_coords2, seeg_systems, NN_flag2] = nifti_values(mni_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'); 

    end
    
end