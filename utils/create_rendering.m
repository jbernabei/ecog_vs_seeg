function [] = create_rendering(rendering_filename, mni_coordinates, adj_matrices, patient_names, varargin)
    % This script performs rendering and saves the output file 
    
    % BrainNet viewer is a dependency
    
    %% Rendering code
    % make matrix
    
    num_patients = length(mni_coordinates);
    
    for pt = 1:num_patients
        pt_name = patient_names{pt};
        
        if rem(nargin,4) == 0
            color_data = -1*ones(size(mni_coordinates{pt},1),1);
            size_data = ones(size(mni_coordinates{pt},1),1);
        else
            color_data = varargin{1};
            size_data = varargin{2};
        end
    
        data_struct = [mni_coordinates{pt}, color_data, size_data];
        
        pt_adj = adj_matrices{pt}(1).data

        % save .node file
        dlmwrite('mni_rendering.node',data_struct,'delimiter',' ','precision',5);
        save('mni_rendering.edge','pt_adj','-ascii');
        
        out_path = sprintf('output/patient_specific/%s/',pt_name);

        % define what to save rendering as & path
        save_to = strcat(out_path, rendering_filename);

        % perform rendering
        BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','mni_rendering.node','mni_rendering.edge','ecog_seeg_render.mat',save_to);

        % delete node file
        delete mni_rendering.node mni_rendering.edge
    
    end

end