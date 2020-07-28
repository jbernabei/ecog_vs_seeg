function [all_metric, pt_metric] = test_EZ(adj_matrices,resected_elecs,metric_type)
    strcmp(metric_type,'node_strength')
    strcmp(metric_type,'betweenness_centrality')
    strcmp(metric_type,'control_centrality')
end