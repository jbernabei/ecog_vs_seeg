function [new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, mod_method)

    num_patients = length(adj_matrices);
    
    for pt = 1:num_patients
        if strcmp(mod_method, 'no_WM')
            %% keep all electrodes except WM
            
            remove_inds = [];
            
            for i = 1:length(patient_roi{pt})
                if patient_roi{pt}(i) == 9171 || patient_roi{pt}(i) == 0
                    remove_inds = [remove_inds, i];
                end
            end
            
            % do adj
            new_adj = adj_matrices{pt};
            new_adj(remove_inds, :) = [];
            new_adj(:, remove_inds) = [];
            
            % do coords
            new_coords = mni_coordinates{pt};
            new_coords(remove_inds, :) = [];
            
            % do roi
            
            % do resected elecs
            
            % assign everything into new structures
              
        elseif strcmp(mod_method, 'min_ROI')
            %% get rid of WM and reduce to ROI-level nodes
            
            % find unique brain regions
            unique_roi = unique(patient_roi{pt});
            
            a = 0;
            
            % loop through unique brain regions
            for i = 1:length(unique_roi)
        
                % find which nodes are in these
                nodes_1 = find(strcmp(patient_roi{pt}, unique_roi(i)));

                % find centroid of these regions
                centroid_1 = mean(mni_coordinates{pt}(nodes_1,:),1);
                pt_centroid{pt}.data(i,:) = centroid_1;
               
                % old code
                 % check if its a resected region
                if sum(sum(nodes_1'==resected_elecs{pt}))>0
                    a = a+1;
                    resect_region(s).data(a) = i;
                end

                % loop through regions again
                for j = 1:length(unique_roi)

                    % find second set of nodes
                    nodes_2 = find(strcmp(patient_roi{pt}, unique_roi(j)));

                    % find second set of centroidss
                    centroid_2 = mean(mni_coordinates{pt}(nodes_2,:),1);
                    if i==j
                    else
                        for k = 1:5
                            reduced_adj(s).freq(k).data(i,j) = mean(mean(adj_matrices{pt}(f).data(nodes_1,nodes_2)));
                        end 
                    end
                end
                
            end
            
            % need new adj
            
            % need new coords
            
            % need new roi
            
            % need new res_elecs
            
        end
        
    end
    
end