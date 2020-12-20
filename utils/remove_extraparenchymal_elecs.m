function [patient_struct] = remove_extraparenchymal_elecs(raw_patient_struct)
    
    % loop through patient structure 
    num_pts = length(raw_patient_struct);
    for pt = 1:num_pts
        % extract roi
        roi_list = raw_patient_struct(pt).roi;
        
        unlocalized_regions = (roi_list==0);
        
        if sum(unlocalized_regions)
            
        fprintf('removing %d unlocalized regions\n',sum(unlocalized_regions))    
            
        resect_bool = zeros(1,length(roi_list));
        resect_bool(raw_patient_struct(pt).resect) = 1;
        
        % now we must modify each item, being the connnectivity,
        % coordinates, roi, and resected electrodes to adequately eliminate
        % unlocalized regions.
        raw_patient_struct(pt).roi(find(unlocalized_regions)) = [];
        raw_patient_struct(pt).coords(find(unlocalized_regions),:) = [];
        
        for i = 1:5
            raw_patient_struct(pt).conn(i).data(:,find(unlocalized_regions)) = [];
            raw_patient_struct(pt).conn(i).data(find(unlocalized_regions),:) = [];
        end
        
        resect_bool(find(unlocalized_regions)) = [];
        
        raw_patient_struct(pt).resect = find(resect_bool);
        
        else
        end
    end
    
    patient_struct = raw_patient_struct;
end