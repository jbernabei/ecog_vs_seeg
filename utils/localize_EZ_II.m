function [targets_1, targets_2] = localize_EZ_II(W, num_targets, mni_coords)

    num_nodes = size(mni_coords,1);
    
    node_str = sum(W);

    % loop through nodes
    for i = 1:num_nodes
        this_node_coords = mni_coords(i,:);
        
        % find all dists
        for j = 1:num_nodes
            dist(j)= sqrt(sum((this_node_coords-mni_coords(j,:)).^2));
        end
        
        [B,I] = sort(dist,'ascend');
        
        node_group(i).data = I(1:(num_targets));
        
        nodestr_resect(i) = sum(node_str(node_group(i).data));
        
        nodeconn_resect(i) = mean(mean(W(node_group(i).data,node_group(i).data)));
    end
    
    % find resection w/ greatest total node strength resected
    [y, which_resection_1] = max(nodestr_resect);
    
    targets_1 = node_group(which_resection_1).data;
    
    % find resection w/ greatest total node strength resected
    [y, which_resection_2] = max(nodeconn_resect);
    
    targets_2 = node_group(which_resection_2).data;
    
end