function [curve, all_dist, all_conn,conn_data] = compute_distance_regression(adj_matrices, mni_coordinates, resected_elecs)
    
    all_conn = []; % initialize elec x 5 freq
    all_dist = [];

    % get number of patients
    num_patients = length(adj_matrices);
    
    % loop through patients
    for pt = 1:num_patients
        
        clear pt_conn
        clear pt_dist
        
        % get number of electrodes
        num_elecs = size(mni_coordinates{pt},1);
        
        % assign into mni coordinates
        mni_coords = mni_coordinates{pt};
        
        % double loop through electrodes and calculate distance matrix
        for i = 1:num_elecs
            for j = 1:num_elecs
                % have to add in resecte elecs here
                pt_dist(i,j) = sqrt(sum((mni_coords(i,:)-mni_coords(j,:)).^2));
            end
        end
        
        all_pt_dist = pt_dist(:);
        
        for f = 1:5
            adj = adj_matrices{pt}(f).data;
            pt_conn(:,f) = adj(:);
        end
        
        all_conn = [all_conn;pt_conn];
        all_dist = [all_dist;all_pt_dist];
        
    end
    
    % extract self connections from data
    self_conn = find(all_dist == 0);
    all_dist(self_conn) = [];
    all_conn(self_conn,:) = [];
    
    % do the curve fitting
    for f = 1:5
        curve(f).data = fit(all_dist,all_conn(:,f),'power2');
        %[fitresult, gof] = powerFit(all_dist,all_conn);
        %curve(f).data = gof;
    end

    % we need to go through each patient and perform the correction to adj
    % matrices based on this distance regression.
    for pt = 1:num_patients
        clear pt_conn
        clear pt_dist
        clear expected_conn
        
        % get number of electrodes
        num_elecs = size(mni_coordinates{pt},1);
        
        % assign into mni coordinates
        mni_coords = mni_coordinates{pt};
        
        % double loop through electrodes and calculate distance matrix
        for i = 1:num_elecs
            for j = 1:num_elecs
                % have to add in resecte elecs here
                pt_dist(i,j) = sqrt(sum((mni_coords(i,:)-mni_coords(j,:)).^2));
                for f = 1:5
                    expected_conn(i,j) = (curve(f).data.a.*pt_dist(i,j).^curve(f).data.b)+curve(f).data.c;
                    new_conn(f).data(i,j) = adj_matrices{pt}(f).data(i,j)-expected_conn(i,j);
                    if i==j
                        new_conn(f).data(i,j) = 0;
                    end
                end
            end
        end
        
        for f = 1:5
        % clip so everything is between 0 and 1
%         new_conn(f).data(new_conn(f).data>1) = 1;
%         new_conn(f).data(new_conn(f).data<0) = 0;

            new_conn(f).data(isinf(new_conn(f).data)) = 0;
        end
        
        
        conn_data(pt).conn = new_conn;
    end
    
    
    
end