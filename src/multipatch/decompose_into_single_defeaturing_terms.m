% DECOMPOSE_INTO_SINGLE_DEFEATURING_TERMS: transforms some defeaturing data
% from the notation of [BCV2022] to the data structures used in this code
% to avoid code repetition. 
%
% CALL:
% 
%   defeaturing_data_into_parts = decompose_into_single_defeaturing_terms(defeaturing_data,...
%                                       problem_data, [problem_data_0, problem_data_tildeFp])
%
% INPUT: (see notation in [BCV2022])
%
%     defeaturing_data:               a structure with data of the defeatured problem. It contains the fields:
%        - gamma_n_sides_in_omega:    ids of the boundaries corresponding to gamma_n in the mesh on Omega
%        - gamma_n_sides_in_omega0:   a structure with data on the boundary gamma_n in the mesh on Omega_0. It contains the fields:
%           -- patch:                 ids of the patches in which lie the sides corresponding to the boundary gamma_n
%           -- local_side_in_patch    ids local to the patch gamma_n_sides_in_omega0.patch of the boundaries corresponding to gamma_n
%       AND/OR
%        - gamma_r_sides_in_omega:    ids of the boundaries corresponding to gamma_r in the mesh on Omega
%        - gamma_r_sides_in_tildeFp:  a structure with data on the boundary gamma_r in the mesh on \tilde F_p. It contains the fields:
%           -- patch:                 ids of the patches in which lie the sides corresponding to the boundary gamma_r
%           -- local_side_in_patch:   ids local to the patch gamma_r_sides_in_tildeFp.patch of the boundaries corresponding to gamma_r
%       AND/OR
%        - gamma_0p_sides_in_omega0:  ids of the boundaries corresponding to gamma_0p in the mesh on Omega_0
%        - gamma_0p_sides_in_tildeFp: ids of the boundaries corresponding to gamma_0p in the mesh on \tilde F_p
%     problem_data:                   a structure with data of the exact problem in \Omega. It contains the fields:
%        - g:                         function for the exact Neumann boundary condition 
%        - (optional) c_diff:         diffusion coefficient
%     (optional) problem_data_0:      a structure with data of the defeatured problem in \Omega_0. It contains the fields:
%        - g:                         function for the defeatured Neumann boundary condition 
%     (optional) problem_data_tildeFp: a structure with data of the extension problem in \tilde F_p. It contains the fields:
%        - g:                         function for the Neumann boundary conditions imposed on the Neumann boundaries of /tilde F_p
%
% OUTPUT:
%
%     defeaturing_data_into_parts: a structure with data considering each boundary contribution to the defeaturing error. It contains the fields:
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


function defeaturing_data_into_parts = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data, varargin)

    if nargin > 2
        problem_data_0 = varargin{1};
        if nargin > 3
            problem_data_tildeFp = varargin{2};
        end
    end

    if ~isfield(problem_data, 'c_diff')
        problem_data.c_diff = @(varargin) ones(size(varargin{1}));
    end

    defeaturing_data_into_parts = defeaturing_data;

    if isfield(defeaturing_data, 'gamma_n_sides_in_omega')
        defeaturing_data_negative.c_diff = problem_data.c_diff;
        defeaturing_data_negative.neumann_function = @(varargin) problem_data.g(varargin{1:end-1}, ...
            defeaturing_data.gamma_n_sides_in_omega(varargin{end}));
        defeaturing_data_negative.defeatured_boundaries_in_defeatured_domain = defeaturing_data.gamma_n_sides_in_omega0;
        defeaturing_data_negative.defeatured_boundaries_in_exact_domain = defeaturing_data.gamma_n_sides_in_omega;

        defeaturing_data_into_parts.gamma_n = defeaturing_data_negative;
    end

    if isfield(defeaturing_data, 'gamma_r_sides_in_omega')
        defeaturing_data_extension.c_diff = problem_data.c_diff;
        defeaturing_data_extension.neumann_function = @(varargin) problem_data_tildeFp.g(varargin{1:end-1}, ...
            defeaturing_data.gamma_r_sides_in_omega(varargin{end}));
        defeaturing_data_extension.defeatured_boundaries_in_defeatured_domain = defeaturing_data.gamma_r_sides_in_tildeFp;
        defeaturing_data_extension.defeatured_boundaries_in_exact_domain = defeaturing_data.gamma_r_sides_in_omega;

        defeaturing_data_into_parts.gamma_r = defeaturing_data_extension;
    end

    if isfield(defeaturing_data, 'gamma_0p_sides_in_omega0')
        defeaturing_data_positive.c_diff = problem_data.c_diff;
        defeaturing_data_positive.neumann_function = @(varargin) -problem_data_0.g(varargin{1:end-1}, ...
            defeaturing_data.gamma_0p_sides_in_omega0(varargin{end})); 
        defeaturing_data_positive.defeatured_boundaries_in_defeatured_domain = defeaturing_data.gamma_0p_sides_in_tildeFp;
        defeaturing_data_positive.defeatured_boundaries_in_exact_domain = defeaturing_data.gamma_0p_sides_in_omega0;

        defeaturing_data_into_parts.gamma_0p = defeaturing_data_positive;
    end

end
