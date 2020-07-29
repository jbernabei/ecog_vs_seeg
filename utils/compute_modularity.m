function [pt_modules, pt_q_vals, pt_pc] = compute_modularity(adj_matrices)
   
    % get number of patients
    num_patients = length(adj_matrices);
    
    % loop through patients
    for pt = 1:num_patients
        
        % loop through frequency bands 
        for f = 1:5
            
            % extract adjacency matrix
            patient_adj = adj_matrices{pt}(f).data;
            
            % compute modularityy (scale factor gamma = 1)
            [S,Q] = modularity_und(patient_adj,1);
            
            % assign modules and Q values
            pt_modules{pt}.data(:,f) = S;
            pt_q_vals(pt,f) = Q;
            
            % calculate nodal participation coefficient
            pt_pc{pt}.data(:,f)=participation_coef(patient_adj,S);
        end
        
    end
    
end