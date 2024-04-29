% ESTIMATE_DEFEATURING_ERROR_H1S: computes the estimator of the energy norm
% of the defeaturing error given in [BCV2022, AC2023+] for a trimmed geometry with
% (possibly multiple) negative features.
%
% CALL:
% 
%   [estimators, measures_of_boundary] = estimate_defeaturing_error_H1s(defeaturing_data, space_0, u_0, msh_trimmed)
%
% INPUT:
%
%     defeaturing_data: a structure with data considering the boundary contributions to the defeaturing error. It contains the fields:
%        - defeatured_boundaries_in_exact_domain:  ids of the boundaries corresponding to gamma in the (trimmed) mesh on Omega
%        - c_diff:                                 exact diffusion coefficient
%        - neumann_function:                       exact Neumann function g, for which boundaries are numbered with respect to the number of sides of the considered defeatured boundary
%     space_0:          space object in Omega_0 that defines the discrete space
%     u_0:              vector of dof weights respresenting the defeatured solution u_0 in Omega_0
%     msh_trimmed:      mesh object for Omega that defines the quadrature rule
%
% OUTPUT:
%
%     estimators:             value of the defeaturing error estimator components corresponding to each defeatured boundary
%     measure_of_boundaries:  measure of each of the different defeatured boundaries (called gamma in [BCV2022, AC2023+])
%
% Copyright (C) 2024 Ondine Chanon
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

function [estimators, measure_of_boundaries] = estimate_defeaturing_error_H1s(defeaturing_data, space_0, u_0, msh_trimmed)

    number_of_boundary_sides = length(defeaturing_data.defeatured_boundaries_in_exact_domain); 

    estimators = zeros(number_of_boundary_sides, 1);
    measure_of_boundaries = zeros(number_of_boundary_sides, 1);
    
    for side_id = 1:number_of_boundary_sides
        neumann_error_non_trimmed = 0;
        weights_times_jacdet_non_trimmed = 0;
        neumann_error_trimmed = 0;
        weights_times_jacdet_trimmed = 0;

        norm_of_neumann_error_squared = 0;
        integral_of_neumann_error = 0;

        side = defeaturing_data.defeatured_boundaries_in_exact_domain(side_id);
        bnd_id = msh_trimming_boundary_find(msh_trimmed.reparam, side);
        param_side = msh_trimmed.reparam.boundaries(bnd_id).parametric_side;

        non_trimmed_elems_bd_ids = msh_trimmed.reparam.boundaries(bnd_id).non_trimmed_elem_bd_ids;
        reparam_elems_bd_ids = msh_trimmed.reparam.boundaries(bnd_id).reparam_elem_bd_ids;

        if ~isempty(non_trimmed_elems_bd_ids) 
            local_indices_non_trimmed = get_boundary_indices(param_side, msh_trimmed.msh_cart.nel_dir, non_trimmed_elems_bd_ids);
            msh_non_trimmed = msh_eval_boundary_side (msh_trimmed.mesh_of_level(ilev).msh_cart, param_side, local_indices_non_trimmed);
            msh_non_trimmed_fi = msh_boundary_side_from_interior(msh_trimmed.msh_cart, param_side);
            sp = space_0.space_untrimmed.constructor(msh_non_trimmed_fi);

            msh_non_trimmed_struct = msh_evaluate_element_list(msh_non_trimmed_fi, local_indices_non_trimmed);
            sp_struct = sp_evaluate_element_list (sp, msh_non_trimmed_struct, 'value', false, 'gradient', true);
            % Pass from the untrimmed to the trimmed numbering
            sp_struct.ndof = space_0.ndof;
            sp_struct.connectivity = space_0.global_to_active(sp_struct.connectivity);

            [gradient_u_0_non_trimmed, F_non_trimmed] = sp_eval_msh(u_0, sp_struct, msh_non_trimmed_struct, 'gradient');
            gradient_u_0_dot_normal_non_trimmed = reshape(sum(gradient_u_0_non_trimmed .* msh_non_trimmed.normal, 1), msh_non_trimmed.nqn, msh_non_trimmed.nel);

            quadrature_points_non_trimmed = cell(msh_trimmed.rdim, 1);
            for idim = 1:msh_trimmed.rdim
                quadrature_points_non_trimmed{idim} = reshape(F_non_trimmed(idim,:), msh_non_trimmed.nqn, msh_non_trimmed.nel);
            end
            c_diff = defeaturing_data.c_diff(quadrature_points_non_trimmed{:});
            exact_neumann = defeaturing_data.neumann_function(quadrature_points_non_trimmed{:}, side_id);

            neumann_error_non_trimmed = exact_neumann - c_diff .* gradient_u_0_dot_normal_non_trimmed;
            weights_times_jacdet_non_trimmed = msh_non_trimmed.quad_weights .* msh_non_trimmed.jacdet;

            norm_of_neumann_error_squared = norm_of_neumann_error_squared ...
                + sum(sum(weights_times_jacdet_non_trimmed .* neumann_error_non_trimmed.^2));
            integral_of_neumann_error = integral_of_neumann_error ...
                + sum(sum(weights_times_jacdet_non_trimmed .* neumann_error_non_trimmed));
            measure_of_boundaries(side_id) = measure_of_boundaries(side_id) + sum(sum(weights_times_jacdet_non_trimmed));
        end

        if ~isempty(reparam_elems_bd_ids)
            msh_tiles_bd = msh_evaluate_boundary_tiles(msh_trimmed, side, true);
            sp_tiles = sp_evaluate_tiles(space_0, msh_tiles_bd, 'gradient', true);

            [gradient_u_0_trimmed, F_trimmed] = sp_eval_msh(u_0, sp_tiles, msh_tiles_bd, 'gradient');
            gradient_u_0_dot_normal_trimmed = reshape(sum(gradient_u_0_trimmed .* msh_tiles_bd.normal, 1), [], msh_tiles_bd.nel);

            quadrature_points_trimmed = cell(msh_trimmed.rdim, 1);
            for idim = 1:msh_trimmed.rdim
                quadrature_points_trimmed{idim} = reshape(F_trimmed(idim,:), [], msh_tiles_bd.nel);
            end
            c_diff = defeaturing_data.c_diff(quadrature_points_trimmed{:});
            exact_neumann = defeaturing_data.neumann_function(quadrature_points_trimmed{:}, side_id);

            neumann_error_trimmed = exact_neumann - c_diff .* gradient_u_0_dot_normal_trimmed;
            weights_times_jacdet_trimmed = msh_tiles_bd.quad_weights .* msh_tiles_bd.jacdet;

            norm_of_neumann_error_squared = norm_of_neumann_error_squared ...
                + sum(sum(weights_times_jacdet_trimmed .* neumann_error_trimmed.^2));
            integral_of_neumann_error = integral_of_neumann_error ...
                + sum(sum(weights_times_jacdet_trimmed .* neumann_error_trimmed));
            measure_of_boundaries(side_id) = measure_of_boundaries(side_id) + sum(sum(weights_times_jacdet_trimmed));
        end

        average_of_d_gamma = integral_of_neumann_error ./ measure_of_boundaries(side_id);

        norm_of_neumann_error_meines_average_squared = ...
            sum(sum(weights_times_jacdet_non_trimmed .* (neumann_error_non_trimmed - average_of_d_gamma).^2))...
            + sum(sum(weights_times_jacdet_trimmed .* (neumann_error_trimmed - average_of_d_gamma).^2));
    
        power_defeaturing = 1 / (msh_trimmed.rdim - 1);
        power_compatibility = msh_trimmed.rdim / (msh_trimmed.rdim - 1);
        if msh_trimmed.rdim == 2
            c_gamma_squared = max(abs(log(measure_of_boundaries(side_id))), lambertw(1)); 
        else
            c_gamma_squared = 1;
        end
    
        estimators(side_id) = sqrt(measure_of_boundaries(side_id).^power_defeaturing .* norm_of_neumann_error_meines_average_squared...
                         + c_gamma_squared * measure_of_boundaries(side_id).^power_compatibility .* average_of_d_gamma.^2);

    end

end
