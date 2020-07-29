function [non_res_metric, res_metric, pt_metric] = test_EZ(adj_matrices,resected_elecs,metric_type)
    
    % calculate number of patients
    num_patients = length(adj_matrices);
    
    % loop through patients
    for pt = 1:num_patients
        
        % extract resected electrode indices
        pt_res_elecs = resected_elecs{pt};
        
        % loop through frequency bands
        for f = 1:5
            
            % extract adjacency matrix
            patient_adj = adj_matrices{pt}(f).data;

            % select which metric we are using
            if strcmp(metric_type,'node_strength')
                pt_metric{pt}.data(:,f) = zscore(sum(patient_adj));
            elseif strcmp(metric_type,'betweenness_centrality')
                pt_metric{pt}.data(:,f) = zscore(betweenness_wei(patient_adj));
            elseif strcmp(metric_type,'control_centrality')
                pt_metric{pt}.data(:,f) = zscore(control_centrality(patient_adj));
            end
           
        end
        
        % calculate mean of the metric across resected electrodes
        res_metric(pt,:) = mean(pt_metric{pt}.data(pt_res_elecs,:),1);

        % extract non-resected electrode data
        non_res = pt_metric{pt}.data;
        non_res(pt_res_elecs,:) = [];

        % calculate mean across non-resected electrodes
        non_res_metric(pt,:) = mean(non_res);
            
    end

end