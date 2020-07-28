function [connectivity_all, connectivity_pt] = analyze_connectivity(adj_matrices)
    num_patients = length(adj_matrices);
    
    connectivity_all = [];
    
    for pt = 1:num_patients
        for f = 1:5
        patient_adj = adj_matrices{pt}(f).data;
        connectivity_pt{pt}.data(:,f) = patient_adj(:);
        end
    end
    
    connectivity_all = [connectivity_all; connectivity_pt{pt}.data];
    
    connectivity_all(find(connectivity_all==0)) = [];
end