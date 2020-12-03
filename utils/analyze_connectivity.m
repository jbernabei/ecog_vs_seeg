function [connectivity_all, connectivity_pt] = analyze_connectivity(adj_matrices,patient_roi)
    num_patients = length(adj_matrices);
    
    connectivity_all = [];
    
    for pt = 1:num_patients
        pt
        for f = 1:5
        patient_adj = adj_matrices{pt}(f).data;
        
        % get patient's WM inds
        patient_WM_inds = find(patient_roi{pt}==9171);
        
        size(patient_adj)
        
        patient_adj(patient_WM_inds,:) = [];
        patient_adj(:,patient_WM_inds) = [];
        
        connectivity_pt{pt}.data(:,f) = patient_adj(:);
        end
        
        connectivity_all = [connectivity_all; connectivity_pt{pt}.data];
    end
    
    
    
    %connectivity_all(find(connectivity_all==0)) = [];
end