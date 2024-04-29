% CREATE_LOOPS_FROM_HOLES: Given the centers and radiuses of circular or rectangular 
% holes in a rectangle geometry, create the loops parameter to execute the
% trimming process with the function ref_trimming_reparameterization_2D.
%
% USAGE:
%
%   loops = create_loops_from_holes(centers, radiuses, types, trimming_option, labels, domain_bounds)
%
% INPUT:
%
%  centers:         centers of the holes, matrix of dimension 2 x #holes
%  radiuses:        radiuses of the holes, vector of length #holes if
%           types==CIRCLE or matrix of dimension 2 x #holes if types==RECTANGLE 
%           (in the latter case, radius = length of the edge / 2).
%  types:           either CIRCLE (= 0), or RECTANGLE (= 1), vector of length #holes
%  trimming_option: either 'param' (cut operation in the parametric domain) or 'physical' (Euclidean cut operation)
%  labels:          label of each trimming curve, vector of length #holes
%  domain_bounds:   bounds of the rectangular domain, [x_min, x_max; y_min, y_max]
%
% OUTPUT:
%
%  loops:    the trimming loops corresponding to the holes, see ref_trimming_reparameterization_2D
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

function loops = create_loops_from_holes(centers, radiuses, types, trimming_option, labels, domain_bounds)

% Types of holes
CIRCLE = 0;
RECTANGLE = 1;

if nargin < 6
    domain_bounds = [0 1; 0 1];
end

if strcmp(trimming_option, 'physical')
    
    loops = cell(size(radiuses, 1), 1);
    for ii = 1:size(radiuses, 1)
        % define the closed loop
        if types(ii) == CIRCLE
            loops{ii}.segments{1}.curve = nrbcirc(radiuses(ii), centers(ii, :));
            loops{ii}.segments{1}.label = labels(ii);
            loops{ii}.tool_type = 1;
        elseif types(ii) == RECTANGLE
            loops{ii}.segments{1}.curve = nrbline([centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)],...
                [centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)]);
            loops{ii}.segments{1}.label = labels(ii);
            loops{ii}.segments{2}.curve = nrbline([centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)],...
                [centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)]);
            loops{ii}.segments{2}.label = labels(ii);
            loops{ii}.segments{3}.curve = nrbline([centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)],...
                [centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)]);
            loops{ii}.segments{3}.label = labels(ii);
            loops{ii}.segments{4}.curve = nrbline([centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)],...
                [centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)]);
            loops{ii}.segments{4}.label = labels(ii);
            loops{ii}.tool_type = 1;
        end
    end
    
elseif strcmp(trimming_option, 'param')
    
    feat_by_bd = cell(5,1);
    for ifeat = 1:size(centers, 1)
        dim2 = 1 + types(ifeat);
        if centers(ifeat, 1) < domain_bounds(1,1) + radiuses(ifeat, 1)
            feat_by_bd{1} = [feat_by_bd{1}, ifeat];
        elseif centers(ifeat, 1) > domain_bounds(1,2) - radiuses(ifeat, 1)
            feat_by_bd{2} = [feat_by_bd{2}, ifeat];
        elseif centers(ifeat, 2) < domain_bounds(2,1) + radiuses(ifeat, dim2)
            feat_by_bd{3} = [feat_by_bd{3}, ifeat];
        elseif centers(ifeat, 2) > domain_bounds(2,2) - radiuses(ifeat, dim2)
            feat_by_bd{4} = [feat_by_bd{4}, ifeat];
        else % internal features
            feat_by_bd{5} = [feat_by_bd{5}, ifeat];
        end
    end
    
    loops = cell(length(feat_by_bd{5})+1, 1);
    
    % Internal features
    if ~isempty(feat_by_bd{5})
        for ifeat = 1:length(feat_by_bd{5})
            ii = feat_by_bd{5}(ifeat);
            if types(ii) == CIRCLE
                loops{ifeat}.segments{1}.curve = nrbcirc(radiuses(ii), centers(ii, :));
                loops{ifeat}.segments{1}.label = labels(ii);
            else % types(ii) == RECTANGLE
                loops{ifeat}.segments{1}.curve = nrbline([centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)],...
                    [centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)]);
                loops{ifeat}.segments{1}.label = labels(ii);
                loops{ifeat}.segments{2}.curve = nrbline([centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)],...
                    [centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)]);
                loops{ifeat}.segments{2}.label = labels(ii);
                loops{ifeat}.segments{3}.curve = nrbline([centers(ii, 1)+radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)],...
                    [centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)]);
                loops{ifeat}.segments{3}.label = labels(ii);
                loops{ifeat}.segments{4}.curve = nrbline([centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)+radiuses(ii, 2)],...
                    [centers(ii, 1)-radiuses(ii, 1), centers(ii, 2)-radiuses(ii, 2)]);
                loops{ifeat}.segments{4}.label = labels(ii);
            end
            loops{ifeat}.tool_type = 3; 
        end
    end
    
    icount = 1;
    loopid = length(feat_by_bd{5}) + 1;
    
    % Features on boundary 1
    if ~isempty(feat_by_bd{1})
        [~, id] = sort(centers(feat_by_bd{1}, 2));
        id = feat_by_bd{1}(id);

        y1 = zeros(length(id), 1);
        y2 = y1;
        theta = y1;
        for ii = 1:length(id)
            if types(id(ii)) == CIRCLE
                theta(ii) = acos( (domain_bounds(1,1)-centers(id(ii), 1)) / radiuses(id(ii)) );
                y2(ii) = radiuses(id(ii)) * sin(theta(ii)) + centers(id(ii), 2);
                y1(ii) = 2*centers(id(ii), 2) - y2(ii);
            else
                y1(ii) = centers(id(ii), 2)-radiuses(id(ii), 2);
                y2(ii) = centers(id(ii), 2)+radiuses(id(ii), 2);
            end
        end

        loops{loopid}.segments{icount}.curve = nrbline(domain_bounds(:,1)', [domain_bounds(1,1), y1(1)]);
        icount = icount + 1;
        for ii = 1:length(id)
            if types(id(ii)) == RECTANGLE
                loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,1), y1(ii)], ...
                    [centers(id(ii), 1)+radiuses(id(ii), 1), y1(ii)]);
                loops{loopid}.segments{icount}.label = labels(id(ii));
                loops{loopid}.segments{icount+1}.curve = nrbline([centers(id(ii), 1)+radiuses(id(ii), 1), y1(ii)], ...
                                       [centers(id(ii), 1)+radiuses(id(ii), 1), y2(ii)]);
                loops{loopid}.segments{icount+1}.label = labels(id(ii));
                loops{loopid}.segments{icount+2}.curve = nrbline([centers(id(ii), 1)+radiuses(id(ii), 1), y2(ii)], ...
                    [domain_bounds(1,1), y2(ii)]);
                loops{loopid}.segments{icount+2}.label = labels(id(ii));
                
                icount = icount + 4;

            elseif types(id(ii)) == CIRCLE
                loops{loopid}.segments{icount}.curve = nrbcirc(radiuses(id(ii)), centers(id(ii), :), -theta(ii), theta(ii));
                loops{loopid}.segments{icount}.label = labels(id(ii));
                
                icount = icount + 2;
            end

            if ii ~= length(id)
                loops{loopid}.segments{icount-1}.curve = nrbline([domain_bounds(1,1), y2(ii)], [domain_bounds(1,1), y1(ii+1)]);
            else
                loops{loopid}.segments{icount-1}.curve = nrbline([domain_bounds(1,1), y2(ii)], [domain_bounds(1,1), domain_bounds(2,2)]);
            end
            loops{loopid}.segments{icount-1}.label = 1;
        end
    else
        loops{loopid}.segments{icount}.curve = nrbline(domain_bounds(:,1)', [domain_bounds(1,1), domain_bounds(2,2)]);
        loops{loopid}.segments{icount}.label = 1; 
        icount = icount + 1;
    end
    
    % Features on boundary 4
    if ~isempty(feat_by_bd{4})
        [~, id] = sort(centers(feat_by_bd{4}, 1));
        id = feat_by_bd{4}(id);

        x1 = zeros(length(id), 1);
        x2 = x1;
        theta = x1;
        for ii = 1:length(id)
            if types(id(ii)) == CIRCLE
                theta(ii) = asin( (domain_bounds(2,2)-centers(id(ii), 2)) / radiuses(id(ii)) );
                x2(ii) = radiuses(id(ii)) * cos(theta(ii)) + centers(id(ii), 1);
                x1(ii) = 2*centers(id(ii), 1) - x2(ii);
            else
                x1(ii) = centers(id(ii), 1)-radiuses(id(ii), 1);
                x2(ii) = centers(id(ii), 1)+radiuses(id(ii), 1);
            end
        end

        loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,1), domain_bounds(2,2)], [x1(1), domain_bounds(2,2)]);
        icount = icount + 1;
        for ii = 1:length(id)
            if types(id(ii)) == RECTANGLE
                loops{loopid}.segments{icount}.curve = nrbline([x1(ii), domain_bounds(2,2)],...
                    [x1(ii), centers(id(ii), 2)-radiuses(id(ii), 2)]);
                loops{loopid}.segments{icount}.label = labels(id(ii));
                loops{loopid}.segments{icount+1}.curve = nrbline([x1(ii), centers(id(ii), 2)-radiuses(id(ii), 2)], ...
                                       [x2(ii), centers(id(ii), 2)-radiuses(id(ii), 2)]);
                loops{loopid}.segments{icount+1}.label = labels(id(ii));
                loops{loopid}.segments{icount+2}.curve = nrbline([x2(ii), centers(id(ii), 2)-radiuses(id(ii), 2)], ...
                    [x2(ii), domain_bounds(2,2)]);
                loops{loopid}.segments{icount+2}.label = labels(id(ii));
                
                icount = icount + 4;

            elseif types(id(ii)) == CIRCLE
                loops{loopid}.segments{icount}.curve = nrbcirc(radiuses(id(ii)), centers(id(ii), :), -pi-theta(ii), theta(ii));
                loops{loopid}.segments{icount}.label = labels(id(ii));
                
                icount = icount + 2;
            end

            if ii~=length(id)
                loops{loopid}.segments{icount-1}.curve = nrbline([x2(ii), domain_bounds(2,2)], [x1(ii+1), domain_bounds(2,2)]);
            else
                loops{loopid}.segments{icount-1}.curve = nrbline([x2(ii), domain_bounds(2,2)], domain_bounds(:,2)');
            end
            loops{loopid}.segments{icount-1}.label = 4;
        end
    else
        loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,1), domain_bounds(2,2)], domain_bounds(:,2)');
        loops{loopid}.segments{icount}.label = 4; 
        icount = icount + 1;
    end
    
    % Features on boundary 2
    if ~isempty(feat_by_bd{2})
        [~, id] = sort(centers(feat_by_bd{2}, 2));
        id = feat_by_bd{2}(id(end:-1:1));

        y1 = zeros(length(id), 1);
        y2 = y1;
        theta = y1;
        for ii = 1:length(id)
            if types(id(ii)) == CIRCLE
                theta(ii) = acos( (domain_bounds(1,2)-centers(id(ii), 1)) / radiuses(id(ii)) );
                y1(ii) = radiuses(id(ii)) * sin(theta(ii)) + centers(id(ii), 2);
                y2(ii) = 2*centers(id(ii), 2) - y1(ii);
            else
                y1(ii) = centers(id(ii), 2)+radiuses(id(ii), 2);
                y2(ii) = centers(id(ii), 2)-radiuses(id(ii), 2);
            end
        end

        loops{loopid}.segments{icount}.curve = nrbline(domain_bounds(:,2)', [domain_bounds(1,2), y1(1)]);
        icount = icount + 1;
        for ii = 1:length(id)
            if types(id(ii)) == RECTANGLE
                loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,2), y1(ii)],...
                    [centers(id(ii), 1)-radiuses(id(ii), 1), y1(ii)]);
                loops{loopid}.segments{icount}.label = labels(id(ii));
                loops{loopid}.segments{icount+1}.curve = nrbline([centers(id(ii), 1)-radiuses(id(ii), 1), y1(ii)], ...
                                       [centers(id(ii), 1)-radiuses(id(ii), 1), y2(ii)]);
                loops{loopid}.segments{icount+1}.label = labels(id(ii));
                loops{loopid}.segments{icount+2}.curve = nrbline([centers(id(ii), 1)-radiuses(id(ii), 1), y2(ii)],...
                    [domain_bounds(1,2), y2(ii)]);
                loops{loopid}.segments{icount+2}.label = labels(id(ii));
                
                icount = icount + 4;
            elseif types(id(ii)) == CIRCLE
                loops{loopid}.segments{icount}.curve = nrbcirc(radiuses(id(ii)), centers(id(ii), :), theta(ii), 2*pi-theta(ii));
                loops{loopid}.segments{icount}.label = labels(id(ii));
            
                icount = icount + 2;
            end

            if ii~=length(id)
                loops{loopid}.segments{icount-1}.curve = nrbline([domain_bounds(1,2), y2(ii)], [domain_bounds(1,2), y1(ii+1)]);
            else
                loops{loopid}.segments{icount-1}.curve = nrbline([domain_bounds(1,2), y2(ii)], [domain_bounds(1,2), domain_bounds(2,1)]);
            end
            loops{loopid}.segments{icount-1}.label = 2;
        end
    else
        loops{loopid}.segments{icount}.curve = nrbline(domain_bounds(:,2)', [domain_bounds(1,2), domain_bounds(2,1)]);
        loops{loopid}.segments{icount}.label = 2; 
        icount = icount + 1;
    end
    
    % Features on boundary 3
    if ~isempty(feat_by_bd{3})
        [~, id] = sort(centers(feat_by_bd{3}, 1));
        id = feat_by_bd{3}(id(end:-1:1));

        x1 = zeros(length(id), 1);
        x2 = x1;
        theta = x1;
        for ii = 1:length(id)
            if types(id(ii)) == CIRCLE
                theta(ii) = asin( (domain_bounds(2,1)-centers(id(ii), 2)) / radiuses(id(ii)) );
                x1(ii) = radiuses(id(ii)) * cos(theta(ii)) + centers(id(ii), 1);
                x2(ii) = 2*centers(id(ii), 1) - x1(ii);
            else
                x1(ii) = centers(id(ii), 1)+radiuses(id(ii), 1);
                x2(ii) = centers(id(ii), 1)-radiuses(id(ii), 1);
            end
        end

        loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,2), domain_bounds(2,1)], [x1(1), domain_bounds(2,1)]);
        icount = icount + 1;
        for ii = 1:length(id)
            if types(id(ii)) == RECTANGLE
                loops{loopid}.segments{icount}.curve = nrbline([x1(ii), domain_bounds(2,1)],...
                    [x1(ii), centers(id(ii), 2)+radiuses(id(ii), 2)]);
                loops{loopid}.segments{icount}.label = labels(id(ii));
                loops{loopid}.segments{icount+1}.curve = nrbline([x1(ii), centers(id(ii), 2)+radiuses(id(ii), 2)], ...
                                       [x2(ii), centers(id(ii), 2)+radiuses(id(ii), 2)]);
                loops{loopid}.segments{icount+1}.label = labels(id(ii));
                loops{loopid}.segments{icount+2}.curve = nrbline([x2(ii), centers(id(ii), 2)+radiuses(id(ii), 2)], ...
                    [x2(ii), domain_bounds(2,1)]);
                loops{loopid}.segments{icount+2}.label = labels(id(ii));
                
                icount = icount + 4;
            elseif types(id(ii)) == CIRCLE
                loops{loopid}.segments{icount}.curve = nrbcirc(radiuses(id(ii)), centers(id(ii), :), theta(ii), pi-theta(ii));
                loops{loopid}.segments{icount}.label = labels(id(ii));

                icount = icount + 2;
            end

            if ii~=length(id)
                loops{loopid}.segments{icount-1}.curve = nrbline([x2(ii), domain_bounds(2,1)], [x1(ii+1), domain_bounds(2,1)]);
            else
                loops{loopid}.segments{icount-1}.curve = nrbline([x2(ii), domain_bounds(2,1)], domain_bounds(:,1)');
            end
            loops{loopid}.segments{icount-1}.label = 3;
        end
    else
        loops{loopid}.segments{icount}.curve = nrbline([domain_bounds(1,2), domain_bounds(2,1)], domain_bounds(:,1)');
        loops{loopid}.segments{icount}.label = 3; 
    end
    
    loops{loopid}.tool_type = 3;
    
end

end

