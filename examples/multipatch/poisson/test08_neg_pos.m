clear
close all
clc

% -- Test reproducing [AC2023+], Section 6.1.2, geometry \Omega_\delta, \delta > 0 --
% -- For delta = 0, see test07_complex.m ; for delta < 0, see test06_complex.m --

%% Set problem data
% -- Choose below the desired delta > 0 to test --
% delta = 2e-1;
delta = 2e-4; 
% --

feature_side_length = 0.1;

filename = 'results/test08_neg_pos';
saveIt = false;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 40 * cos(pi * x) + 10 * cos(5 * pi * x);

[problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data);
defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_Fp);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32]; 
method_data.nquad = [5 5];


%% Main

% 1) BUILD GEOMETRY
[srf_0, srf, srf_Fp] = build_geometry(feature_side_length, delta); 
problem_data.geo_name = srf;
problem_data_0.geo_name = srf_0;
problem_data_Fp.geo_name = srf_Fp;

% 2) SOLVE THE EXACT PROBLEM
[omega, msh, space, u] = mp_solve_laplace_generalized(problem_data, method_data);

% 3a) SOLVE THE DEFEATURED PROBLEM
[omega_0, msh_0, space_0, u_0] = mp_solve_laplace_generalized(problem_data_0, method_data);
problem_data_Fp.h = extract_boundary_dofs(msh_0, space_0, u_0, defeaturing_data.gamma_0p_sides_in_omega0);

% 3b) SOLVE THE EXTENSION PROBLEM
[Fp, msh_Fp, space_Fp, u_0tilde] = mp_solve_laplace_generalized(problem_data_Fp, method_data);

% 4a) COMPUTE THE DEFEATURING ESTIMATOR
[estimator, measure_of_gamma_n, measure_of_gamma_0p] = ...
    estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0, msh_Fp, space_Fp, u_0tilde);

% 4b) COMPUTE THE DEFEATURING ERROR
[error_H1s, error_H1s_Omega_star, error_H1s_Fp] = ...
    defeaturing_error_H1s(msh_0, space_0, u_0, defeaturing_data.omega_star_patches_in_omega0, ...
                          msh, space, u, defeaturing_data.omega_star_patches_in_omega, ...
                          msh_Fp, space_Fp, u_0tilde, defeaturing_data.Fp_patches_in_tildeFp);
error_H1s_boundary_representation = ...
    defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, ...
                          msh_0, space_0, u_0, msh_Fp, space_Fp, u_0tilde);


%% Display and save the results
if saveIt
    save(filename, 'feature_side_length', 'delta', 'error_H1s', 'error_H1s_Omega_star', 'error_H1s_Fp',...
        'error_H1s_boundary_representation', 'measure_of_gamma_n', 'measure_of_gamma_0p', 'estimator')
end

fprintf('For delta = %e, \n', delta)
fprintf('    * Estimator                           = %e \n', estimator)
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_H1s)
fprintf('    * Effectivity index                   = %f \n', estimator / error_H1s)


%% Auxiliary functions
function [srf_0, srf, srf_Fpos] = build_geometry(epsilon, delta)
    extension_factor = 4; 

    srf_tot(1) = nrb4surf([0, 0], [0.5 - epsilon - delta / 2, 0], ...
                          [0, 1 - epsilon], [0.5 - epsilon - delta / 2, 1 - epsilon]);
    srf_tot(1) = nrbkntins(srf_tot(1), {(1 - extension_factor * epsilon) / (1 - epsilon), ...
                                         1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(2) = nrb4surf([0.5 - epsilon - delta / 2, 0], [0.5 - delta / 2, 0], ...
                          [0.5 - epsilon - delta / 2, 1 - epsilon], [0.5 - delta / 2, 1 - epsilon]);
    srf_tot(2) = nrbkntins(srf_tot(2), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(3) = nrb4surf([0.5 - delta / 2, 0], [0.5 + delta / 2, 0], ...
                          [0.5 - delta / 2, 1 - epsilon], [0.5 + delta / 2, 1 - epsilon]);
    srf_tot(3) = nrbkntins(srf_tot(3), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(4) = nrb4surf([0.5 + delta / 2, 0], [0.5 + epsilon + delta / 2, 0], ...
                          [0.5 + delta / 2, 1 - epsilon], [0.5 + epsilon + delta / 2, 1 - epsilon]);
    srf_tot(4) = nrbkntins(srf_tot(4), {min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(5) = nrb4surf([0.5 + epsilon + delta / 2, 0], [1, 0], ...
                          [0.5 + epsilon + delta / 2, 1 - epsilon], [1, 1 - epsilon]);
    srf_tot(5) = nrbkntins(srf_tot(5), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    
    srf_tot(6) = nrb4surf([0, 1 - epsilon], [0.5 - epsilon - delta / 2, 1 - epsilon], ...
                          [0, 1], [0.5 - epsilon - delta / 2, 1]);
    srf_tot(6) = nrbkntins(srf_tot(6), {(1 - extension_factor * epsilon) / (1 - epsilon), ...
                                         1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(7) = nrb4surf([0.5 - epsilon - delta / 2, 1 - epsilon], [0.5 - delta / 2, 1 - epsilon], ...
                          [0.5 - epsilon - delta / 2, 1], [0.5 - delta / 2, 1]);
    srf_tot(7) = nrbkntins(srf_tot(7), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(8) = nrb4surf([0.5 - delta / 2, 1 - epsilon], [0.5 + delta / 2, 1 - epsilon], ...
                          [0.5 - delta / 2, 1], [0.5 + delta / 2, 1]);
    srf_tot(8) = nrbkntins(srf_tot(8), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(9) = nrb4surf([0.5 + delta / 2, 1 - epsilon], [0.5 + epsilon + delta / 2, 1 - epsilon], ...
                          [0.5 + delta / 2, 1], [0.5 + epsilon + delta / 2, 1]);
    srf_tot(9) = nrbkntins(srf_tot(9), {min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(10) = nrb4surf([0.5 + epsilon + delta / 2, 1 - epsilon], ...
                           [1, 1 - epsilon], [0.5 + epsilon + delta / 2, 1], [1, 1]);
    srf_tot(10) = nrbkntins(srf_tot(10), {1-(1-extension_factor*epsilon)/(1-epsilon), ...
                                          1-(extension_factor-1)*epsilon/2});
    
    srf_tot(11) = nrb4surf([0.5 - epsilon - delta / 2, 1], [0.5 - delta / 2, 1], ...
                           [0.5 - epsilon - delta / 2, 1 + epsilon], [0.5 - delta / 2, 1 + epsilon]);
    srf_tot(11) = nrbkntins(srf_tot(11), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), 0.5});
    

    srf_0 = srf_tot(1:10);
    srf_Fpos = srf_tot(11);
    srf = srf_tot([1:8 10:11]);
end

function [problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    problem_data_Fp = problem_data;
    problem_data_Fp = rmfield(problem_data_Fp, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 6 7 9:18]; 
    problem_data.drchlt_sides = [2 3 4 5 8];

    % Simplified problem
    problem_data_0.nmnn_sides = [1 6 8:14];
    problem_data_0.drchlt_sides = [2 3 4 5 7];
    
    % Extension problem
    problem_data_Fp.nmnn_sides = [1 2 4]; 
    problem_data_Fp.drchlt_sides = 3; 
end

function defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_Fp)
    defeaturing_data.omega_star_patches_in_omega0 = [1:8 10];
    defeaturing_data.omega_star_patches_in_omega = 1:9;
    defeaturing_data.Fp_patches_in_tildeFp = 1; 

    defeaturing_data.gamma_n_sides_in_omega = [6 11 13]; 
    defeaturing_data.gamma_n_sides_in_omega0.patch = [4 8 10];
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = [4 2 1];

    defeaturing_data.gamma_0p_sides_in_omega0 = 10;
    defeaturing_data.gamma_0p_sides_in_tildeFp = 3; 

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data, problem_data_0, problem_data_Fp);
end
