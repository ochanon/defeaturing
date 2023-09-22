% ESTIMATE_SINGLE_DEFEATURING_ERROR_TERM_H1S: computes one term of the
% estimator of the energy norm of the defeaturing error given in [BCV2022],
% correponding to one of the boundaries gamma_n, gamma_r or gamma_0p. 
%
% CALL:
% 
%   [estimator, measure_of_boundary] = estimate_single_defeaturing_error_term_H1s(defeaturing_data, ...
%                                                       defeatured_msh, defeatured_space, defeatured_u)
%
% INPUT:
%
%     defeaturing_data:   a structure with data considering one boundary contribution to the defeaturing error. It contains the fields:
%         - c_diff:                                      exact diffusion coefficient
%         - neumann_function:                            exact Neumann function g, for which boundaries are numbered with respect to the number of sides of the considered defeatured boundary
%         - defeatured_boundaries_in_defeatured_domain:  a structure with data on the defeatured boundary in the mesh on the defeatured geometry. It contains the fields:
%           -- patch:                 ids of the patches in which lie the sides corresponding to the defeatured boundary
%           -- local_side_in_patch:   ids local to the patch defeatured_boundaries_in_defeatured_domain.patch of the sides corresponding to the defeatured boundary
%     defeatured_msh:     mesh object of the considered defeatured geometry, that defines the quadrature rule
%     defeatured_space:   space object in the considered defeatured geometry, that defines the discrete space
%     defeatured_u:       vector of dof weights respresenting the defeatured solution, in the considered defeatured geometry
%     
% OUTPUT:
%
%     estimator:           value of the term of the defeaturing error estimator corresponding to the considered defeatured boundary
%     measure_of_boundary: measure of the considered defeatured boundary
%
% Copyright (C) 2023 Ondine Chanon
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.


function [estimator, measure_of_boundary] = estimate_single_defeaturing_error_term_H1s(defeaturing_data, ...
                                                defeatured_msh, defeatured_space, defeatured_u)
    norm_of_neumann_error_squared = 0;
    integral_of_neumann_error = 0;
    measure_of_boundary = 0;

    if ~isstruct(defeaturing_data.defeatured_boundaries_in_defeatured_domain)
        number_of_boundary_sides = length(defeaturing_data.defeatured_boundaries_in_defeatured_domain); 
        boundary_sides_in_defeatured_domain = zeros(number_of_boundary_sides, 1);
        patches_of_boundary_sides_in_defeatured_domain = zeros(number_of_boundary_sides, 1);
        side_id = 0;
        for side = defeaturing_data.defeatured_boundaries_in_defeatured_domain
            for bnd_side = 1:defeatured_msh.boundaries(side).nsides
                side_id = side_id + 1;
                patches_of_boundary_sides_in_defeatured_domain(side_id) = defeatured_msh.boundaries(side).patches(bnd_side);
                boundary_sides_in_defeatured_domain(side_id) = defeatured_msh.boundaries(side).faces(bnd_side);
            end
        end
    else
        patches_of_boundary_sides_in_defeatured_domain = defeaturing_data.defeatured_boundaries_in_defeatured_domain.patch;
        boundary_sides_in_defeatured_domain = defeaturing_data.defeatured_boundaries_in_defeatured_domain.local_side_in_patch;
    end

    number_of_boundary_sides = length(boundary_sides_in_defeatured_domain); 
    neumann_error = cell(number_of_boundary_sides, 1);
    weights_times_jacdet = cell(number_of_boundary_sides, 1);
    
    for side_id = 1:number_of_boundary_sides
        patch = patches_of_boundary_sides_in_defeatured_domain(side_id); 
        local_side = boundary_sides_in_defeatured_domain(side_id);

        [gradient_u_0, msh_0_side] = sp_eval_msh_on_patch_side(defeatured_msh, defeatured_space, defeatured_u, patch, local_side, 'gradient');
        gradient_u_0_dot_normal = reshape(sum (gradient_u_0 .* msh_0_side.normal, 1), msh_0_side.nqn, msh_0_side.nel);

        quadrature_points = cell(defeatured_msh.rdim, 1);
        for idim = 1:defeatured_msh.rdim
            quadrature_points{idim} = reshape(msh_0_side.geo_map(idim, :, :), msh_0_side.nqn, msh_0_side.nel);
        end
        c_diff = defeaturing_data.c_diff(quadrature_points{:});
        exact_neumann = defeaturing_data.neumann_function(quadrature_points{:}, side_id);

        neumann_error{side_id} = exact_neumann - c_diff .* gradient_u_0_dot_normal;
        weights_times_jacdet{side_id} = msh_0_side.quad_weights .* msh_0_side.jacdet;

        norm_of_neumann_error_squared = norm_of_neumann_error_squared ...
            + sum(sum(weights_times_jacdet{side_id} .* neumann_error{side_id}.^2));
        integral_of_neumann_error = integral_of_neumann_error ...
            + sum(sum(weights_times_jacdet{side_id} .* neumann_error{side_id}));
        measure_of_boundary = measure_of_boundary + sum(sum(weights_times_jacdet{side_id}));
    end

    average_of_d_gamma = integral_of_neumann_error ./ measure_of_boundary;

    norm_of_neumann_error_meines_average_squared = 0;
    for side_id = 1:number_of_boundary_sides
        norm_of_neumann_error_meines_average_squared = norm_of_neumann_error_meines_average_squared ...
            + sum(sum(weights_times_jacdet{side_id} .* (neumann_error{side_id} - average_of_d_gamma).^2));
    end

    power_defeaturing = 1 / (defeatured_msh.rdim - 1);
    power_compatibility = defeatured_msh.rdim / (defeatured_msh.rdim - 1);
    if defeatured_msh.rdim == 2
        c_gamma_squared = max(abs(log(measure_of_boundary)), lambertw(1)); 
    else
        c_gamma_squared = 1;
    end

    estimator = sqrt(measure_of_boundary.^power_defeaturing .* norm_of_neumann_error_meines_average_squared...
                     + c_gamma_squared * measure_of_boundary.^power_compatibility .* average_of_d_gamma.^2);
end
