function [distances_all, distances_pt] = compute_interelectrode_distances(mni_coordinates)
    num_patients = length(mni_coordinates);
    
    distances_all = [];
    
    for pt = 1:num_patients
        mni_coords = mni_coordinates{pt};
        
        pt_num_elecs = size(mni_coords,1);
        
        distances_pt{pt} = NaN*ones(pt_num_elecs);
        
        for i = 1:pt_num_elecs
            for j = 1:pt_num_elecs
                if i==j
                    
                else
                    elec1 = mni_coords(i,:);
                    elec2 = mni_coords(j,:);
                    distances_pt{pt}(i,j) = sqrt(sum((elec1-elec2).^2));
                end
            end
        end
        
        distances_all = [distances_all; distances_pt{pt}(:)];
        
    end
    
    distances_all = rmmissing(distances_all);
    
end