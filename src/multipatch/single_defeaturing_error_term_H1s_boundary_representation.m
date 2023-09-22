% SINGLE_DEFEATURING_ERROR_TERM_H1s_BOUNDARY_REPRESENTATION: Evaluate one term of the defeaturing error in 
% the H1-semi norm computed from the boundary representation formula, i.e.,
%   ||\nabla (u-u_0)||_L^2(\Omega) = int_\gamma (g - \nabla u_0 \dot n)(u - u_0) ds,
% with \gamma being either \gamma_n, \gamma_r or \gamma_0p. 
% The notation used corresponds to the one in [BCV2022].
%
% CALL:
% 
%   error_H1s = single_defeaturing_error_term_H1s_boundary_representation(defeaturing_data, exact_msh, exact_space, exact_u, ...
%                                                                                defeatured_msh, defeatured_space, defeatured_u)
%
% INPUT:
%
%     defeaturing_data:   a structure with data considering one boundary contribution to the defeaturing error. It contains the fields:
%         - c_diff:                                      exact diffusion coefficient
%         - neumann_function:                            exact Neumann function g, for which boundaries are numbered with respect to the number of sides of the considered defeatured boundary
%         - defeatured_boundaries_in_exact_domain:       ids of the sides corresponding to the defeatured boundary in the mesh on Omega
%         - defeatured_boundaries_in_defeatured_domain:  a structure with data on the defeatured boundary in the mesh on the defeatured geometry. It contains the fields:
%           -- patch:                 ids of the patches in which lie the sides corresponding to the defeatured boundary
%           -- local_side_in_patch:   ids local to the patch defeatured_boundaries_in_defeatured_domain.patch of the sides corresponding to the defeatured boundary
%     exact_msh:          mesh object of the considered exact geometry, that defines the quadrature rule
%     exact_space:        space object in the considered exact geometry, that defines the discrete space
%     exact_u:            vector of dof weights respresenting the exact solution, in the considered exact geometry
%     defeatured_msh:     mesh object of the considered defeatured geometry, that defines the quadrature rule
%     defeatured_space:   space object in the considered defeatured geometry, that defines the discrete space
%     defeatured_u:       vector of dof weights respresenting the defeatured solution, in the considered defeatured geometry
%
% OUTPUT:
%
%     error_H1s: value of the term of the defeaturing error in H^1-semi norm, corresponding to the considered defeatured boundary
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


function error_H1s = single_defeaturing_error_term_H1s_boundary_representation(defeaturing_data, exact_msh, exact_space, exact_u, ...
                                                                               defeatured_msh, defeatured_space, defeatured_u)
    if ~isfield(defeaturing_data, 'c_diff')
        defeaturing_data.c_diff = @(varargin) ones(size(varargin{1}));
    end

    number_of_boundary_sides = length(defeaturing_data.defeatured_boundaries_in_exact_domain); 
    boundary_sides_in_exact_domain = zeros(number_of_boundary_sides, 1);
    patches_of_boundary_sides_in_exact_domain = zeros(number_of_boundary_sides, 1);
    side_id = 0;
    for side = defeaturing_data.defeatured_boundaries_in_exact_domain
        for bnd_side = 1:exact_msh.boundaries(side).nsides
            side_id = side_id + 1;
            patches_of_boundary_sides_in_exact_domain(side_id) = exact_msh.boundaries(side).patches(bnd_side);
            boundary_sides_in_exact_domain(side_id) = exact_msh.boundaries(side).faces(bnd_side);
        end
    end

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

    error_H1s_squared = 0;
    for side_id = 1:number_of_boundary_sides
        patch = patches_of_boundary_sides_in_exact_domain(side_id);
        local_side = boundary_sides_in_exact_domain(side_id); 
        value_u = sp_eval_msh_on_patch_side(exact_msh, exact_space, exact_u, patch, local_side, 'value');
        
        patch_0 = patches_of_boundary_sides_in_defeatured_domain(side_id);
        local_side_0 = boundary_sides_in_defeatured_domain(side_id);
        [u_0_evaluation, msh_0_side] = sp_eval_msh_on_patch_side(defeatured_msh, defeatured_space, defeatured_u, patch_0, local_side_0, {'value', 'gradient'});

        value_u_0 = u_0_evaluation{1};
        gradient_u_0 = u_0_evaluation{2};
        gradient_u_0_dot_normal = reshape(sum (gradient_u_0 .* msh_0_side.normal, 1), msh_0_side.nqn, msh_0_side.nel);
    
        quadrature_points = cell(defeatured_msh.rdim, 1);
        for idim = 1:defeatured_msh.rdim
            quadrature_points{idim} = reshape(msh_0_side.geo_map(idim, :, :), msh_0_side.nqn, msh_0_side.nel);
        end
        c_diff = defeaturing_data.c_diff(quadrature_points{:});
        exact_neumann = defeaturing_data.neumann_function(quadrature_points{:}, side_id);

        neumann_error = exact_neumann - c_diff .* gradient_u_0_dot_normal;
        defeaturing_error = value_u - value_u_0;
        weights_times_jacdet = msh_0_side.quad_weights .* msh_0_side.jacdet;

        error_H1s_squared = error_H1s_squared + sum(sum(weights_times_jacdet .* neumann_error .* defeaturing_error));
    end

    error_H1s = sqrt(error_H1s_squared);
end
