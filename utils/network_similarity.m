function [pt_purity, all_purity] = network_similarity(network_modules, mni_coordinates, comparison_type)
    
    % get number of patients
    num_patients = length(mni_coordinates);
    
    for pt = 1:num_patients
        
        clear purity_index
        
        % extract coordinates
        mni_coords = mni_coordinates{pt};
        
        % compute systems
        [~, seven_systems, ~] = nifti_values(mni_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii'); 

        % need to generate spatial indicator vectors for both modules and
        % systems
        for s = 1:7
            ind_vec_sys(pt).data(s,:) = double([seven_systems==s]);
            
            for f = 1:5
                        
            % modules
            pt_modules = network_modules{pt}.data;
        
            
            % get it for correct frequency
            this_freq_module = pt_modules(:,f);
            
                % loop through modules
                for m = 1:length(unique(this_freq_module))
                    ind_vec_mod(pt).patient(f).data(m,:) = double([this_freq_module==m]);
                    
                    JI_val = 1-pdist([ind_vec_sys(pt).data(s,:); ind_vec_mod(pt).patient(f).data(m,:)],'jaccard');
                    
                    % calculating jaccard index
                    JI(pt).data(s,m,f) = JI_val;
                end

            end
            
        end
        
        
        % calculate Jaccard index
        %JI = 1-pdist([x],'jaccard');
        
        % must compare against null model
        
        % each module should have 7 overlap scores (one for each system)
        
        % compute entropy over elements as purity (-sum((s_i)log2(si)))
        % this should have dimensions as rows-> modules, columns ->
        % frequency
        for f = 1:5
            this_freq_module = pt_modules(:,f);
            for m = 1:length(unique(this_freq_module))
                p = squeeze(JI(pt).data(:,m,f)./sum(JI(pt).data(:,m,f)));
                purity_index(m,f) = max(p);
            end
        end
      
        pt_purity{pt} = purity_index;
        
        all_purity(pt,:) = mean(purity_index);
        
    end
    
end