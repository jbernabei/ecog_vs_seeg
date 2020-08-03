function [color_map] = make_colormap(colorstart, colormid, colorend)
    
    for i=1:3
        first_half(:,i) = linspace(colorstart(i),colormid(i),50);
        second_half(:,i) = linspace(colormid(i),colorend(i),50);
    end
    
    color_map = [first_half;second_half];

end