% ESTIMATE_DEFEATURING_ERROR_H1S: computes the estimator of the energy norm
% of the defeaturing error given in [BCV2022]. 
%
% CALL:
% 
%   [estimator, measure_of_gamma_n, measure_of_gamma_r, measure_of_gamma_0p] = ...
%           estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0, [msh_tildeFp, space_tildeFp, u_0tilde])
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
%     msh_0:              mesh object for Omega_0 that defines the quadrature rule
%     space_0:            space object in Omega_0 that defines the discrete space
%     u_0:                vector of dof weights respresenting the defeatured solution u_0 in Omega_0
%     (optional but should be present all three of them, if present:)
%     msh_tildeFp:        mesh object for \tilde F_p that defines the quadrature rule
%     space_tildeFp:      space object in \tilde F_p that defines the discrete space
%     u_0tilde:           vector of dof weights representing the extended defeatured solution \tilde u_0 in \tilde F_p
%     
% OUTPUT:
%
%     estimator:           value of the defeaturing error estimator
%     measure_of_gamma_n:  (output present if gamma_n is not empty) measure of the boundary gamma_n
%     measure_of_gamma_r:  (output present if gamma_r is not empty) measure of the boundary gamma_r
%     measure_of_gamma_0p: (output present if gamma_0p is not empty) measure of the boundary gamma_0p
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


function [estimator, varargout] = estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0, varargin)

    if nargin > 4
        msh_tildeFp = varargin{1};
        space_tildeFp = varargin{2}; 
        u_0tilde = varargin{3};
    end

    out_id = 1; 
    if isfield(defeaturing_data, 'gamma_n')
        [estimator_gamma_n, measure_of_gamma_n] = ...
            estimate_single_defeaturing_error_term_H1s(defeaturing_data.gamma_n, msh_0, space_0, u_0);
        varargout{out_id} = measure_of_gamma_n;
        out_id = out_id + 1;
    else
        estimator_gamma_n = 0;
    end

    if isfield(defeaturing_data, 'gamma_r')
        [estimator_gamma_r, measure_of_gamma_r] = ...
            estimate_single_defeaturing_error_term_H1s(defeaturing_data.gamma_r, msh_tildeFp, space_tildeFp, u_0tilde);
        varargout{out_id} = measure_of_gamma_r;
        out_id = out_id + 1;
    else
        estimator_gamma_r = 0;
    end

    if isfield(defeaturing_data, 'gamma_0p')
        [estimator_gamma_0p, measure_of_gamma_0p] = ...
            estimate_single_defeaturing_error_term_H1s(defeaturing_data.gamma_0p, msh_tildeFp, space_tildeFp, u_0tilde);
        varargout{out_id} = measure_of_gamma_0p;
    else
        estimator_gamma_0p = 0;
    end
    
    estimator = sqrt(estimator_gamma_n.^2 + estimator_gamma_r.^2 + estimator_gamma_0p.^2);
end
