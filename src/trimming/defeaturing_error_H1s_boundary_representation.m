% DEFEATURING_ERROR_H1S_BOUNDARY_REPRESENTATION: Evaluate the defeaturing error in the H1-semi norm,
% for a trimmed geometry with negative features, using the boundary representation of the error in [BCV2022, AC2023+]. 
% The notation used corresponds to the one in [BCV2022, AC2023+].
%
% CALL:
% 
%   error_H1s = defeaturing_error_H1s_boundary_representation(defeaturing_data, space_0, u_0, ...
%                                                             msh_trimmed, sp_trimmed, u)
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
%     space_trimmed:    space object in Omega that defines the discrete space
%     u:                vector of dof weights representing the exact solution u in Omega
%
% OUTPUT:
%
%     error_H1s:        defeaturing error in the H^1 seminorm of Omega
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

function error_H1s = defeaturing_error_H1s_boundary_representation(defeaturing_data, space_0, u_0, ...
                                                                   msh_trimmed, sp_trimmed, u)
    number_of_boundary_sides = length(defeaturing_data.defeatured_boundaries_in_exact_domain); 
    error_H1s_squared = 0;
    
    for side_id = 1:number_of_boundary_sides
        side_in_exact_domain = defeaturing_data.defeatured_boundaries_in_exact_domain(side_id);
        reparam_boundary_index = msh_trimming_boundary_find(msh_trimmed.reparam, side_in_exact_domain);
        parametric_side = msh_trimmed.reparam.boundaries(reparam_boundary_index).parametric_side;

        non_trimmed_elems_bd_ids = msh_trimmed.reparam.boundaries(reparam_boundary_index).non_trimmed_elem_bd_ids;
        reparam_elems_bd_ids = msh_trimmed.reparam.boundaries(reparam_boundary_index).reparam_elem_bd_ids;

        if ~isempty(non_trimmed_elems_bd_ids) 
            local_indices_non_trimmed = get_boundary_indices(parametric_side, msh_trimmed.msh_cart.nel_dir, non_trimmed_elems_bd_ids);
            msh_non_trimmed = msh_eval_boundary_side (msh_trimmed.msh_cart, parametric_side, local_indices_non_trimmed);
            msh_non_trimmed_fi = msh_boundary_side_from_interior(msh_trimmed.msh_cart, parametric_side);

            sp_0 = space_0.space_untrimmed.constructor(msh_non_trimmed_fi);
            sp = sp_trimmed.space_untrimmed.constructor(msh_non_trimmed_fi);

            msh_non_trimmed_struct = msh_evaluate_element_list(msh_non_trimmed_fi, local_indices_non_trimmed);

            sp_0_struct = sp_evaluate_element_list (sp_0, msh_non_trimmed_struct, 'value', true, 'gradient', true);
            sp_0_struct.ndof = space_0.ndof;
            sp_0_struct.connectivity = space_0.global_to_active(sp_0_struct.connectivity);

            sp_struct = sp_evaluate_element_list (sp, msh_non_trimmed_struct, 'value', true);
            sp_struct.ndof = sp_trimmed.ndof;
            sp_struct.connectivity = sp_trimmed.global_to_active(sp_struct.connectivity);

            [u_0_evaluation, F_non_trimmed] = sp_eval_msh(u_0, sp_0_struct, msh_non_trimmed_struct, {'value', 'gradient'});
            value_u_0_non_trimmed = u_0_evaluation{1};
            gradient_u_0_non_trimmed = u_0_evaluation{2};
            gradient_u_0_dot_normal_non_trimmed = reshape(sum(gradient_u_0_non_trimmed .* msh_non_trimmed.normal, 1), msh_non_trimmed.nqn, msh_non_trimmed.nel);
            
            value_u_ex_non_trimmed = sp_eval_msh(u, sp_struct, msh_non_trimmed_struct, 'value');
            
            quadrature_points_non_trimmed = cell(msh_trimmed.rdim, 1);
            for idim = 1:msh_trimmed.rdim
                quadrature_points_non_trimmed{idim} = reshape(F_non_trimmed(idim,:), msh_non_trimmed.nqn, msh_non_trimmed.nel);
            end
            c_diff = defeaturing_data.c_diff(quadrature_points_non_trimmed{:});
            exact_neumann = defeaturing_data.neumann_function(quadrature_points_non_trimmed{:}, side_id);

            neumann_error = exact_neumann - c_diff .* gradient_u_0_dot_normal_non_trimmed;
            defeaturing_error = value_u_ex_non_trimmed - value_u_0_non_trimmed;
            weights_times_jacdet = msh_non_trimmed.quad_weights .* msh_non_trimmed.jacdet;
    
            error_H1s_squared = error_H1s_squared + sum(sum(weights_times_jacdet .* neumann_error .* defeaturing_error));
        end

        if ~isempty(reparam_elems_bd_ids)
            msh_tiles_bd = msh_evaluate_boundary_tiles(msh_trimmed, side_in_exact_domain, true);
            sp_0_tiles = sp_evaluate_tiles(space_0, msh_tiles_bd, 'value', true, 'gradient', true);
            sp_tiles = sp_evaluate_tiles(sp_trimmed, msh_tiles_bd, 'value', true);

            [u_0_evaluation, F_trimmed] = sp_eval_msh(u_0, sp_0_tiles, msh_tiles_bd, {'value', 'gradient'});
            value_u_0_trimmed = u_0_evaluation{1};
            gradient_u_0_trimmed = u_0_evaluation{2};
            gradient_u_0_dot_normal_trimmed = reshape(sum(gradient_u_0_trimmed .* msh_tiles_bd.normal, 1), [], msh_tiles_bd.nel);

            value_u_ex_trimmed = sp_eval_msh(u, sp_tiles, msh_tiles_bd, 'value');

            quadrature_points_trimmed = cell(msh_trimmed.rdim, 1);
            for idim = 1:msh_trimmed.rdim
                quadrature_points_trimmed{idim} = reshape(F_trimmed(idim,:), [], msh_tiles_bd.nel);
            end
            c_diff = defeaturing_data.c_diff(quadrature_points_trimmed{:});
            exact_neumann = defeaturing_data.neumann_function(quadrature_points_trimmed{:}, side_id);

            neumann_error = exact_neumann - c_diff .* gradient_u_0_dot_normal_trimmed;
            defeaturing_error = value_u_ex_trimmed - value_u_0_trimmed;
            weights_times_jacdet = msh_tiles_bd.quad_weights .* msh_tiles_bd.jacdet;
    
            error_H1s_squared = error_H1s_squared + sum(sum(weights_times_jacdet .* neumann_error .* defeaturing_error));
        end
    end

    error_H1s = sqrt(error_H1s_squared);
end
