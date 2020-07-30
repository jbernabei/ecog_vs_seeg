function [new_data] = modify_networks(adj_matrices, mni_coordinates, patient_roi, resected_elecs, mod_method)

    num_patients = length(adj_matrices);
    
    for pt = 1:num_patients
        if strcmp(mod_method, 'no_WM')
            %% keep all electrodes except WM
            % set up blank indices to remove
            remove_inds = [];
            
            % transfer ROI into placeholder
            this_pt_roi = patient_roi{pt};
            
            for i = 1:length(this_pt_roi)
                if this_pt_roi(i) == 9171 || this_pt_roi(i) == 0
                    remove_inds = [remove_inds, i];
                end
            end
            
            remove_inds;
            
            % do adj
            for f = 1:5
                new_adj(f).data = adj_matrices{pt}(f).data;
                num_elecs = size(new_adj(f).data,1); % get number of electrodes
                new_adj(f).data(remove_inds, :) = [];
                new_adj(f).data(:, remove_inds) = [];
            end
            
            % do coords
            new_coords = mni_coordinates{pt};
            new_coords(remove_inds, :) = [];
            
            % do roi
            new_roi = this_pt_roi;
            new_roi(remove_inds) = [];
            
            % do resected elecs
            resected_elec_bool = zeros(num_elecs,1);
            resected_elec_bool(resected_elecs{pt}) = 1;
            resected_elec_bool(remove_inds) = [];
            new_res_elecs = find(resected_elec_bool);
            
            % assign everything into new structures
            new_data(pt).conn = new_adj;
            new_data(pt).coords = new_coords;
            new_data(pt).resect = new_res_elecs;
            new_data(pt).roi = new_roi;
              
        elseif strcmp(mod_method, 'min_ROI')
            %% get rid of WM and reduce to ROI-level nodes
            clear reduced_adj
            pt;
            
            % find unique brain regions
            unique_roi = unique(patient_roi{pt});
            
            a = 0;
            
            % loop through unique brain regions
            for i = 1:length(unique_roi)
        
                % find which nodes are in these
                nodes_1 = find(patient_roi{pt}==unique_roi(i));
                
                if size(resected_elecs{pt},2)==1
                    resected_elecs{pt} = resected_elecs{pt}';
                end
                
                % find centroid of these regions
                centroid_1 = mean(mni_coordinates{pt}(nodes_1,:),1);
                pt_centroid{pt}.data(i,:) = centroid_1;
               
                % old code
                 % check if its a resected region
                if sum(sum(nodes_1'==resected_elecs{pt}))>0
                    a = a+1;
                    resect_region(pt).data(a) = i;
                end

                % loop through regions again
                for j = 1:length(unique_roi)

                    % find second set of nodes
                    nodes_2 = find(strcmp(patient_roi{pt}, unique_roi(j)));

                    % find second set of centroidss
                    centroid_2 = mean(mni_coordinates{pt}(nodes_2,:),1);
                    if i==j
                    else
                        for f = 1:5
                            reduced_adj(f).data(i,j) = mean(mean(adj_matrices{pt}(f).data(nodes_1,nodes_2)));
                        end 
                    end
                end
                
            end
            
            % need new adj
            new_data(pt).conn = reduced_adj;
            
            % need new coords
            new_data(pt).coords = pt_centroid{pt}.data;
            
            % need new roi
            new_data(pt).roi = unique_roi;
            
            % need new res_elecs
            new_data(pt).resect = resect_region(pt).data;
            
        end
        
    end
    
end