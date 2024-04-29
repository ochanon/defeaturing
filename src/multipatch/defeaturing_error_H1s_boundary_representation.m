% DEFEATURING_ERROR_H1S_BOUNDARY_REPRESENTATION: Evaluate the defeaturing error in the H1-semi norm,
% for a multipatch geometry with a single (negative, positive or complex) feature, using the boundary representation 
% of the error in [BCV2022]. The notation used corresponds to the one in [BCV2022].
%
% CALL:
% 
%   error_H1s = defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, ...
%                                                             msh_0, space_0, u_0, [msh_tildeFp, space_tildeFp, u_0tilde])
%
% INPUT:
%
%     defeaturing_data:   a structure with data considering each boundary contribution to the defeaturing error. It contains the fields:
%        - gamma_n:                a structure containing the information regarding the boundary \gamma_n. It contains the fields:
%           -- c_diff:                                      exact diffusion coefficient
%           -- neumann_function:                            exact Neumann function g, for which boundaries are numbered with respect to the number of boundaries in gamma_n
%           -- defeatured_boundaries_in_defeatured_domain:  a structure with data on the boundary gamma_n in the mesh on Omega_0, see defeaturing_data.gamma_n_sides_in_omega0 above.
%           -- defeatured_boundaries_in_exact_domain:       ids of the boundaries corresponding to gamma_n in the mesh on Omega
%       AND/OR
%        - gamma_r:                a structure containing the information
%        regarding the boundary \gamma_r. It contains the fields:
%           -- c_diff:                                      exact diffusion coefficient
%           -- neumann_function:                            exact Neumann function g, for which boundaries are numbered with respect to the number of boundaries in gamma_r
%           -- defeatured_boundaries_in_defeatured_domain:  a structure with data on the boundary gamma_r in the mesh on \tilde F_p, see defeaturing_data.gamma_r_sides_in_tildeFp above.
%           -- defeatured_boundaries_in_exact_domain:       ids of the boundaries corresponding to gamma_r in the mesh on Omega
%       AND/OR
%        - gamma_0p:               a structure containing the information
%        regarding the boundary \gamma_0p. It contains the fields:
%           -- c_diff:                                      exact diffusion coefficient
%           -- neumann_function:                            defeatured Neumann function -g_0, for which boundaries are numbered with respect to the number of boundaries in gamma_0p
%           -- defeatured_boundaries_in_defeatured_domain:  a structure with data on the boundary gamma_0p in the mesh on \tilde F_p, see defeaturing_data.gamma_0p_sides_in_tildeFp above.
%           -- defeatured_boundaries_in_exact_domain:       ids of the boundaries corresponding to gamma_0p in the mesh on Omega_0
%     msh:                          mesh object for Omega that defines the quadrature rule
%     space:                        space object in Omega that defines the discrete space
%     u:                            vector of dof weights representing the exact solution u in Omega
%     msh_0:                        mesh object for Omega_0 that defines the quadrature rule
%     space_0:                      space object in Omega_0 that defines the discrete space
%     u_0:                          vector of dof weights respresenting the defeatured solution u_0 in Omega_0
%     (optional but all three of them should be present, if present:)
%     msh_tildeFp:                  mesh object for \tilde F_p that defines the quadrature rule
%     space_tildeFp:                space object in \tilde F_p that defines the discrete space
%     u_0tilde:                     vector of dof weights representing the extended defeatured solution \tilde u_0 in \tilde F_p
%
% OUTPUT:
%
%     error_H1s:             defeaturing error in the H^1 seminorm of Omega
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


function error_H1s = defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, ...
                                                                   msh_0, space_0, u_0, varargin)
    if nargin > 7
        msh_tildeFp = varargin{1};
        space_tildeFp = varargin{2}; 
        u_0tilde = varargin{3};
    end

    error_H1s_squared = 0;
    if isfield(defeaturing_data, 'gamma_n')
        error_H1s_squared = error_H1s_squared + ...
            single_defeaturing_error_term_H1s_boundary_representation(defeaturing_data.gamma_n,...
                msh, space, u, msh_0, space_0, u_0)^2;
    end
    if isfield(defeaturing_data, 'gamma_r')
        error_H1s_squared = error_H1s_squared + ...
            single_defeaturing_error_term_H1s_boundary_representation(defeaturing_data.gamma_r,...
                msh, space, u, msh_tildeFp, space_tildeFp, u_0tilde)^2;
    end
    if isfield(defeaturing_data, 'gamma_0p')
        u_in_omega0 = restrict_to_omega0(u, space, space_0, defeaturing_data);

        error_H1s_squared = error_H1s_squared + ...
            single_defeaturing_error_term_H1s_boundary_representation(defeaturing_data.gamma_0p,...
                msh_0, space_0, u_in_omega0, msh_tildeFp, space_tildeFp, u_0tilde)^2;
    end

    error_H1s = sqrt(error_H1s_squared);
end

function u_in_omega0 = restrict_to_omega0(u, space, space_0, defeaturing_data)
    omega_star_patches_in_omega = defeaturing_data.omega_star_patches_in_omega;

    u_in_omega0 = zeros(space_0.ndof, 1);
    for i_patch = 1:length(omega_star_patches_in_omega)
        u_in_omega0(space_0.gnum{i_patch}) = u(space.gnum{omega_star_patches_in_omega(i_patch)});
    end
end
