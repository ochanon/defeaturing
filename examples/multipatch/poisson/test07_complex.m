clear
close all
clc

%% Set problem data
% -- To reproduce [BCV2022], Section 6.2.1, geometry \Omega_\varepsilon^5 --
epsilon_values = 1e-2 ./ 2.^(0:6);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 10 * cos(3 * pi * x) .* sin(5 * pi * y); 
problem_data.h = @(x, y, ind) zeros(size(x));
% --

% -- To reproduce [AC2023+], Section 6.1.2, geometry \Omega_\delta with \delta = 0 --
% -- For delta > 0, see test08_neg_pos.m ; for delta < 0, see test06_complex.m --
% eps_values = 0.1; % This value correspond to the length side of the squares.
% problem_data.c_diff = @(x, y) ones(size(x));
% problem_data.g = @(x, y, ind) zeros(size(x));
% problem_data.f = @(x, y) zeros(size(x)); 
% problem_data.h = @(x, y, ind) 40 * cos(pi * x) + 10 * cos(5 * pi * x);
% --

filename = 'results/test07_complex';
saveIt = false;
plotIt = true;

[problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data);
defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_Fp);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32]; 
method_data.nquad = [5 5];


%% Main
number_of_epsilons = numel(epsilon_values);
error_H1s = zeros(1, number_of_epsilons);
error_H1s_Omega_star = zeros(1, number_of_epsilons);
error_H1s_Fp = zeros(1, number_of_epsilons);
error_H1s_boundary_representation = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
measure_of_gamma_n = zeros(1, number_of_epsilons);
measure_of_gamma_0p = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_H1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = epsilon_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_Fp] = build_geometry(epsilon); 
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
    [estimator(iter), measure_of_gamma_n(iter), measure_of_gamma_0p(iter)] = ...
        estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0, msh_Fp, space_Fp, u_0tilde);

    % 4b) COMPUTE THE DEFEATURING ERROR
    [error_H1s(iter), error_H1s_Omega_star(iter), error_H1s_Fp(iter)] = ...
        defeaturing_error_H1s(msh_0, space_0, u_0, defeaturing_data.omega_star_patches_in_omega0, ...
                              msh, space, u, defeaturing_data.omega_star_patches_in_omega, ...
                              msh_Fp, space_Fp, u_0tilde, defeaturing_data.Fp_patches_in_tildeFp);
    error_H1s_boundary_representation(iter) = ...
        defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, ...
                              msh_0, space_0, u_0, msh_Fp, space_Fp, u_0tilde);

    norm_of_u(iter) = error_H1s_in_patches(msh, space, u, 1:msh.npatch, ...
                                           msh, space, zeros(size(u)));
    relative_error_H1s(iter) = error_H1s(iter) / norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'epsilon_values', 'error_H1s', 'error_H1s_Omega_star', 'error_H1s_Fp', ...
        'error_H1s_boundary_representation', 'estimator', 'measure_of_gamma_n', ...
        'measure_of_gamma_0p', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', ...
           epsilon_values, epsilon_values, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf("For epsilon (feature's side length) = %e, \n", epsilon)
fprintf('    * Estimator                           = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_H1s(end))
fprintf('    * Effectivity index                   = %f \n', estimator(end) / error_H1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_F_positive] = build_geometry(epsilon)
    extension_factor = 4; 

    srf_tot(1) = nrb4surf([0, 0], [0.5 - epsilon, 0], [0, 1 - epsilon], [0.5 - epsilon, 1 - epsilon]);
    srf_tot(1) = nrbkntins(srf_tot(1), {[1 - 2 * 0.1 / 9, (1 - extension_factor * epsilon) / (1 - epsilon)], ...
                                        [1 - (extension_factor - 1) * epsilon / 2, 0.9 + 8 * 0.1 / 9]});
    
    srf_tot(2) = nrb4surf([0.5 - epsilon, 0], [0.5, 0], [0.5 - epsilon, 1 - epsilon], [0.5, 1 - epsilon]);
    srf_tot(2) = nrbkntins(srf_tot(2), {[0.1, 0.9], [1 - (extension_factor - 1) * epsilon / 2, 0.9 + 8 * 0.1 / 9]});
    
    srf_tot(3) = nrb4surf([0.5, 0], [0.5 + epsilon, 0], [0.5, 1 - epsilon], [0.5 + epsilon, 1 - epsilon]);
    srf_tot(3) = nrbkntins(srf_tot(3), {[0.1, 0.9], [1 - (extension_factor - 1) * epsilon / 2, 0.9 + 8 * 0.1 / 9]});
    
    srf_tot(4) = nrb4surf([0.5 + epsilon, 0], [1, 0], [0.5 + epsilon, 1 - epsilon], [1, 1 - epsilon]);
    srf_tot(4) = nrbkntins(srf_tot(4), {[2 * 0.1 / 9, 1 - (1 - extension_factor * epsilon) / (1 - epsilon)], ...
                                        [1 - (extension_factor - 1) * epsilon / 2, 0.9 + 8 * 0.1 / 9]});
    
    srf_tot(5) = nrb4surf([0, 1 - epsilon], [0.5 - epsilon, 1 - epsilon], [0, 1], [0.5 - epsilon, 1]);
    srf_tot(5) = nrbkntins(srf_tot(5), {[1 - 2 * 0.1 / 9, (1 - extension_factor * epsilon) / (1 - epsilon)], ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(6) = nrb4surf([0.5 - epsilon, 1 - epsilon], [0.5, 1 - epsilon], [0.5 - epsilon, 1], [0.5, 1]);
    srf_tot(6) = nrbkntins(srf_tot(6), {[0.1, 0.9], 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(7) = nrb4surf([0.5, 1 - epsilon], [0.5 + epsilon, 1 - epsilon], [0.5, 1], [0.5 + epsilon, 1]);
    srf_tot(7) = nrbkntins(srf_tot(7), {[0.1, 0.9], 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(8) = nrb4surf([0.5 + epsilon, 1 - epsilon], [1, 1 - epsilon], [0.5 + epsilon, 1], [1, 1]);
    srf_tot(8) = nrbkntins(srf_tot(8), {[2 * 0.1 / 9, 1 - (1 - extension_factor * epsilon) / (1 - epsilon)], ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(9) = nrb4surf([0.5 - epsilon, 1], [0.5, 1], [0.5 - epsilon, 1 + epsilon], [0.5, 1 + epsilon]);
    srf_tot(9) = nrbkntins(srf_tot(9), {[0.1, 0.9], 0.5});
    

    srf_0 = srf_tot(1:8);
    srf_F_positive = srf_tot(9);
    srf = srf_tot([1:6 8 9]);
end

function [problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    problem_data_Fp = problem_data;
    problem_data_Fp = rmfield(problem_data_Fp, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 5 6 8:16]; 
    problem_data.drchlt_sides = [2 3 4 7];
    % problem_data.gamma_sides = [5 10 11]; % for the negative feature
    % problem_data.omega0_patches = 1:7;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 5 7:12];
    problem_data_0.drchlt_sides = [2 3 4 6];
    % problem_data_0.gamma0_sides = 9; % for the Dirichlet extension
    % problem_data_0.omega_patches = [1:6 8];
    % problem_data_0.gamma_sides = cell(8, 1); % relative to each patch
    % problem_data_0.gamma_sides([3 6 8]) = {4, 2, 1};
    
    % Extension problem
    problem_data_Fp.nmnn_sides = [1 2 4]; 
    problem_data_Fp.drchlt_sides = 3; 
    % problem_data_Fp.gamma0_sides = 3; 
    % problem_data_Fp.gammae_sides = {[]}; % relative to each patch
    % problem_data_Fp.F_patches = 1;
end

function defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_Fp)
    defeaturing_data.omega_star_patches_in_omega0 = [1:6 8];
    defeaturing_data.omega_star_patches_in_omega = 1:7;
    defeaturing_data.Fp_patches_in_tildeFp = 1;

    defeaturing_data.gamma_n_sides_in_omega = [5 10 11];
    defeaturing_data.gamma_n_sides_in_omega0.patch = [3 6 8];
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = [4 2 1];

    defeaturing_data.gamma_0p_sides_in_omega0 = 9; 
    defeaturing_data.gamma_0p_sides_in_tildeFp = 3; 

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data, problem_data_0, problem_data_Fp);
end
