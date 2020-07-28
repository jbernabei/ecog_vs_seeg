function [pt_modules, pt_q_vals, pt_pc] = compute_modularity(adj_matrices)
    % using modularity_und as the method
    num_patients = length(patient_conn);
    
    for pt = 1:num_patients
        for f = 1:5
            patient_adj = adj_matrices{pt}(f).data;
            [S,Q] = modularity_und(patient_adj,1);
            pt_modules{pt}.data(:,f) = S;
            pt_q_vals(pt,f) = Q;
            pt_pc{pt}.data(:,f)=participation_coef(patient_adj,S);
        end
        
    end
    
end