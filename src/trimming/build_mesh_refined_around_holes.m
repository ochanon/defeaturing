% BUILD_MESH_REFINED_AROUND_HOLES: Given the centers and radiuses of holes in the
% square geometry (0, 1)^2, a tensor-product mesh on (0, 1)^2 is built such
% that each one of the holes is in a region with fine elements.
%
% USAGE:
%
%  geo_name = build_mesh_refined_around_holes(num_elements, centers, radiuses, degree)
%
% INPUT:
%
%  num_elements: indicative number of elements on which to build the mesh.
%       More precisely, each hole is embedded in a submesh of 2 * num_elements per
%       side; a refinement zone (called extension) around each hole which is 4 times larger 
%       than the hole is also embedded in a submesh of 2 * num_elements in
%       the direction of the extension; then num_elements^2 cover the
%       remaining space of the geometry (0, 1)^2.
%  centers:      coordinates of the holes' centers, matrix of dimension 2 x #holes
%  radiuses:     radiuses of the holes, vector of length #holes
%  degree:       required global mesh degree (kept the same in both space directions)
%
% OUTPUT:
%
%  geo_name:    
%
% Copyright (C) 2024 Ondine Chanon
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function geo_name = build_mesh_refined_around_holes(num_elements, centers, radiuses, degree)

num_elements_in_hole = 2 * num_elements;
num_elements_in_extension = 2 * num_elements;
extension_factor = 4; 

if length(num_elements) == 1
    num_elements = num_elements * ones(1, size(centers, 2));
    num_elements_in_extension = num_elements_in_extension * ones(1, size(centers, 2));
    num_elements_in_hole = num_elements_in_hole * ones(1, size(centers, 2));
end

knots = cell(1, size(centers, 2));

for idim = 1:size(centers, 2)
    [centers_x, Isort] = unique(centers(:, idim)); % centers sorted
    radiuses_x = radiuses(Isort);
    h_out = 1. / num_elements(idim);

    h_extension = (extension_factor - 1) * radiuses_x ./ num_elements_in_extension(idim); 
    h_hole = radiuses_x ./ num_elements_in_hole(idim); 
    hole_boundaries = [centers_x - radiuses_x; centers_x + radiuses_x];
    half_h_hole_vec = [h_hole / 2 ; h_hole / 2];
    h_extension_vec = [h_extension; h_extension];
    indicators = [ones(length(centers_x), 1); -ones(length(centers_x), 1)];
    
    [hole_boundaries, sort_id] = sort(hole_boundaries);
    indicators = indicators(sort_id);
    half_h_hole_vec = half_h_hole_vec(sort_id);
    h_extension_vec = h_extension_vec(sort_id);
    
    id_first = find(hole_boundaries >= 0, 1, 'first');
    id_last = find(hole_boundaries <= 1, 1, 'last');
    hole_boundaries = hole_boundaries(id_first:id_last);
    indicators = indicators(id_first:id_last);
    half_h_hole_vec = half_h_hole_vec(id_first:id_last);
    h_extension_vec = h_extension_vec(id_first:id_last);
    
    indicators = cumsum(indicators);
    left = id_first - 1;
    
    % get rid of multiple centers
    [hole_boundaries, id_unique] = unique(hole_boundaries);
    half_h_hole_vec_unique = zeros(length(id_unique), 1);
    h_extension_vec_unique = zeros(length(id_unique), 1);
    indicators_unique = zeros(length(id_unique), 1);
    id_unique = [id_unique; length(half_h_hole_vec)];
    for ii = 1:length(id_unique)-1
        if ii ~= length(id_unique)-1
            ids = id_unique(ii):id_unique(ii+1)-1;
        else 
            ids = id_unique(ii):id_unique(ii+1);
        end
        half_h_hole_vec_unique(ii) = max(half_h_hole_vec(ids));
        h_extension_vec_unique(ii) = max(h_extension_vec(ids));
        indicators_unique(ii) = indicators(ids(end));
    end
    indicators = indicators_unique;
    half_h_hole_vec = half_h_hole_vec_unique;
    h_extension_vec = h_extension_vec_unique;
    
    % take care of the boundaries of (0, 1)^2
    if hole_boundaries(1) ~= 0
        indicators = [left; indicators];
        hole_boundaries = [0; hole_boundaries];
        half_h_hole_vec = [0; half_h_hole_vec];
        h_extension_vec = [0; h_extension_vec];
        flag0 = true;
    else
        flag0 = false;
    end
    if indicators(end) == 0 && hole_boundaries(end) ~= 1
        indicators = [indicators; 0];
        hole_boundaries = [hole_boundaries; 1];
        half_h_hole_vec = [half_h_hole_vec; 0];
        h_extension_vec = [h_extension_vec; 0];
        flag1 = true;
    elseif hole_boundaries(end)~=1
        indicators = [indicators; 1];
        hole_boundaries = [hole_boundaries; 1];
        half_h_hole_vec = [half_h_hole_vec; 0];
        h_extension_vec = [h_extension_vec; 0];
        flag1 = true;
    else
        flag1 = false;
    end
    
    % find ids of boundaries that begin a "non hole" zone
    non_hole_ids = find(~indicators(1:end-1));
    
    % do the same for the extension zone
    extension_boundaries = [centers_x - extension_factor * radiuses_x; centers_x + extension_factor * radiuses_x]; 
    indicators_extension = [ones(length(centers_x), 1); -ones(length(centers_x), 1)];
    [extension_boundaries, sortid_extension] = sort(extension_boundaries);
    indicators_extension = indicators_extension(sortid_extension);
    
    id_first = find(extension_boundaries >= 0, 1, 'first');
    id_last = find(extension_boundaries <= 1, 1, 'last');
    extension_boundaries = extension_boundaries(id_first:id_last);
    indicators_extension = indicators_extension(id_first:id_last);
    
    indicators_extension = cumsum(indicators_extension) + id_first - 1;
    left = id_first - 1;
    if min(indicators_extension) < 0
        indicators_extension = indicators_extension - min(indicators_extension);
    end
    
    % get rid of multiple centers
    [extension_boundaries, id_unique_extension] = unique(extension_boundaries);
    indicators_unique_extension = zeros(length(id_unique_extension), 1);
    if ~isempty(id_unique_extension)
        id_unique_extension = [id_unique_extension; id_unique_extension(end)];
        for ii = 1:length(id_unique_extension)-1
            if ii ~= length(id_unique_extension)-1
                ids_extension = id_unique_extension(ii):id_unique_extension(ii+1)-1;
            else
                ids_extension = id_unique_extension(ii):id_unique_extension(ii+1);
            end
            indicators_unique_extension(ii) = indicators_extension(ids_extension(end));
        end
    end
    indicators_extension = indicators_unique_extension;
    
    % take care of boundaries
    if isempty(extension_boundaries) || extension_boundaries(1) ~= 0
        indicators_extension = [left; indicators_extension];
        extension_boundaries = [0; extension_boundaries];
    end
    if ~isempty(indicators_extension) && indicators_extension(end) == 0 && extension_boundaries(end) ~= 1
        indicators_extension = [indicators_extension; 0];
        extension_boundaries = [extension_boundaries; 1];
    elseif extension_boundaries(end) ~= 1
        indicators_extension = [indicators_extension; 1];
        extension_boundaries = [extension_boundaries; 1];
    end
    
    % find ids of boundaries that begin a "non hole" zone
    non_extension_ids = find(~indicators_extension(1:end-1));
    extid_real = [];
    for index = 1:length(non_extension_ids)
        extid_temp = find(hole_boundaries <= extension_boundaries(non_extension_ids(index)));
        extid_real = [extid_real, extid_temp(end)];
    end

    knots_x = [];
    hh_icpr = half_h_hole_vec;
    for icpr = 1:length(hole_boundaries)-1
        if (hole_boundaries(icpr)+half_h_hole_vec(icpr) > hole_boundaries(icpr+1)-half_h_hole_vec(icpr+1))
            half_h_hole_vec([icpr, icpr+1]) = (hole_boundaries(icpr+1) - hole_boundaries(icpr)) / 2;
        end
    end
    
    for icpr = 1:length(hole_boundaries)-1
        if icpr == 1 && ~flag0
            knots_x = [knots_x, 0];
        end
        
        if ~ismember(icpr, non_hole_ids) % inside features
            h_icpr = 2 * min(hh_icpr([icpr, icpr + 1])); 
            if h_icpr == 0 && icpr ~= length(hh_icpr)
                h_icpr = 2 * max(hh_icpr([icpr, icpr + 1])); 
            end
            knots_x = [knots_x, linspace(hole_boundaries(icpr) + half_h_hole_vec(icpr), hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1), ...
                ceil((hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1) - hole_boundaries(icpr) - half_h_hole_vec(icpr)) / h_icpr) + 1)];
            
        elseif ~ismember(icpr, extid_real) % in the extension zone between features or at boundaries
            h_icpr = min(h_extension_vec([icpr, icpr + 1])); 
            if h_icpr == 0
                h_icpr = max(h_extension_vec([icpr, icpr + 1])); 
            end
            knots_x = [knots_x, linspace(hole_boundaries(icpr) + half_h_hole_vec(icpr), hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1), ...
                ceil((hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1) - hole_boundaries(icpr) - half_h_hole_vec(icpr)) / h_icpr) + 1)];
            
        else % extension + outside zones between features or at boundaries
            ids_cprext = find(extension_boundaries > hole_boundaries(icpr) + half_h_hole_vec(icpr) & extension_boundaries < hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1));
            if icpr == 1 && extid_real(1) == 1
                knots_x = [knots_x, linspace(hole_boundaries(icpr) + half_h_hole_vec(icpr), extension_boundaries(ids_cprext(1)), ...
                    ceil((extension_boundaries(ids_cprext(1)) - hole_boundaries(icpr) - half_h_hole_vec(icpr)) / h_out) + 1)];
            else
                h_icpr = h_extension_vec(icpr); 
                if h_icpr == 0 && icpr ~= length(h_extension_vec)
                    h_icpr = h_extension_vec(icpr + 1); 
                elseif h_icpr == 0 && icpr ~= 1
                    h_icpr = h_extension_vec(icpr - 1);
                end
                knots_x = [knots_x, linspace(hole_boundaries(icpr) + half_h_hole_vec(icpr), extension_boundaries(ids_cprext(1)), ...
                    ceil((extension_boundaries(ids_cprext(1)) - hole_boundaries(icpr) - half_h_hole_vec(icpr)) / h_icpr) + 1)];
            end
            knots_x = knots_x(1:end-1);
            
            for icprext = 1:length(ids_cprext)-1
                if ismember(ids_cprext(icprext), non_extension_ids) % outside zone
                    knots_x = [knots_x, linspace(extension_boundaries(ids_cprext(icprext)), extension_boundaries(ids_cprext(icprext) + 1),...
                        ceil((extension_boundaries(ids_cprext(icprext) + 1) - extension_boundaries(ids_cprext(icprext))) / h_out) + 1)]; 
                    knots_x = knots_x(1:end-1);
                    
                else % extension zone
                    h_icpr = h_extension_vec(icpr); 
                    if h_icpr == 0 && icpr ~= length(h_extension_vec)
                        h_icpr = h_extension_vec(icpr + 1); 
                    elseif h_icpr == 0 && icpr ~= 1
                        h_icpr = h_extension_vec(icpr - 1);
                    end
                    knots_x = [knots_x, linspace(extension_boundaries(ids_cprext(icprext)), extension_boundaries(ids_cprext(icprext) + 1),...
                        ceil((extension_boundaries(ids_cprext(icprext) + 1) - extension_boundaries(ids_cprext(icprext))) / h_icpr) + 1)];
                    knots_x = knots_x(1:end-1);
                    
                end
            end
            if ismember(ids_cprext(end), non_extension_ids)
                knots_x = [knots_x, linspace(extension_boundaries(ids_cprext(end)), hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1), ...
                    ceil((hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1) - extension_boundaries(ids_cprext(end))) / h_out) + 1)];
            else
                h_icpr = h_extension_vec(icpr + 1); 
                if h_icpr == 0
                    h_icpr = h_extension_vec(icpr);
                end
                knots_x = [knots_x, linspace(extension_boundaries(ids_cprext(end)), hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1), ...
                    ceil((hole_boundaries(icpr + 1) - half_h_hole_vec(icpr + 1) - extension_boundaries(ids_cprext(end))) / h_icpr) + 1)];
            end
        end
        
        if icpr+1 == length(hole_boundaries) && ~flag1
            knots_x = [knots_x, 1];
        end
    end
    
    knots{idim} = knots_x;
end

knots{1} = knots{1}(2:end-1);
knots{2} = knots{2}(2:end-1);

geo_name = nrbsquare([0, 0], 1, 1);
geo_name = nrbdegelev(geo_name, degree - 1);
geo_name = nrbkntins(geo_name, knots);

end
