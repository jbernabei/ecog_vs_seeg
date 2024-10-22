function [targets_1, dice_scores] = localize_EZ(adj_matrices, mni_coordinates, resected_elecs, metric_type)
    num_patients = length(adj_matrices);
    
    for pt = 1:num_patients

        for f = 1:5
            patient_adj = adj_matrices{pt}(f).data;
            
            if strcmp(metric_type,'node_strength')
                localization_metric = sum(patient_adj);
            elseif strcmp(metric_type,'betweenness_centrality')
                localization_metric = -1*betweenness_wei(patient_adj);
            elseif strcmp(metric_type,'control_centrality')
                localization_metric = control_centrality(patient_adj);
            end
            
            % below is old code
            num_nodes = size(patient_adj,1);
            mni_coords = mni_coordinates{pt};
            % loop through nodes
            for i = 1:num_nodes
                this_node_coords = mni_coords(i,:);

                % find all dists
                for j = 1:num_nodes
                    dist(pt).data(j)= sqrt(sum((this_node_coords-mni_coords(j,:)).^2));
                end

                
                [B,I] = sort(dist(pt).data,'ascend');
                
                num_targets = length(resected_elecs{pt});

                node_group(i).data = I(1:(num_targets));

                metric_resect(i) = sum(localization_metric(node_group(i).data));

                nodeconn_resect(i) = mean(mean(patient_adj(node_group(i).data,node_group(i).data)));
            end

            % find resection w/ greatest total node strength resected
            [y, which_resection] = max(metric_resect);

            targets_1 = node_group(which_resection).data;
            
            %targets_1
            %resected_elecs{pt}
            % compute dice
            dice_scores(pt,f) = dice(targets_1,resected_elecs{pt});
            
        end  
        
    end
end