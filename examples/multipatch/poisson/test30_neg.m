clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.2.2, geometry \Omega_\varepsilon^1 --

%% Set problem data
epsilon_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test30_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.f = @(x, y, z) 10 * cos(3 * pi * x) .* sin(5 * pi * y) .* sin(7 * pi * z); 
problem_data.h = @(x, y, z, ind) zeros(size(x));

[problem_data, problem_data_0] = set_boundary_conditions(problem_data);
defeaturing_data = set_defeaturing_data(problem_data);

method_data.degree = [3 3 3];
method_data.regularity = [2 2 2];
method_data.nsub = [8 8 8];
method_data.nquad = [5 5 5];


%% Main
number_of_epsilons = numel(epsilon_values);
error_H1s = zeros(1, number_of_epsilons);
error_H1s_boundary_representation = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
measure_of_gamma = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_H1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = epsilon_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = build_geometry(epsilon);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_generalized(problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_generalized(problem_data_0, method_data);

    % 4a) COMPUTE THE DEFEATURING ESTIMATOR
    [estimator(iter), measure_of_gamma(iter)] = estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0);

    % 4b) COMPUTE THE DEFEATURING ERROR
    error_H1s(iter) = defeaturing_error_H1s(msh_0, space_0, u_0, defeaturing_data.omega_star_patches_in_omega0, ...
                                            msh, space, u);
    error_H1s_boundary_representation(iter) = ...
        defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, msh_0, space_0, u_0);

    norm_of_u(iter) = error_H1s_in_patches(msh, space, u, 1:msh.npatch, ...
                                           msh, space, zeros(size(u)));
    relative_error_H1s(iter) = error_H1s(iter) / norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'epsilon_values', 'error_H1s', 'error_H1s_boundary_representation', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', epsilon_values, epsilon_values.^(3/2)*0.3, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^{3/2}', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(epsilon)
    extension_factor = 4;

    srf_0(1) = nrb4surf([0, 0], [0.5 - epsilon / 2, 0], [0, 1 - epsilon], [0.5 - epsilon / 2, 1 - epsilon]);
    srf_0(1) = nrbextrude(srf_0(1), [0, 0, epsilon]);
    srf_0(1) = nrbkntins(srf_0(1), {(0.5 - extension_factor * epsilon / 2) / (0.5 - epsilon / 2), ...
                                    (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2), []});
    
    srf_0(2) = nrb4surf([0.5 - epsilon / 2, 0], [0.5 + epsilon / 2, 0], ...
                        [0.5 - epsilon / 2, 1 - epsilon], [0.5 + epsilon / 2, 1 - epsilon]);
    srf_0(2) = nrbextrude(srf_0(2), [0, 0, epsilon]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2), []});
    
    srf_0(3) = nrb4surf([0.5 + epsilon / 2, 0], [1, 0], [0.5 + epsilon / 2, 1 - epsilon], [1, 1 - epsilon]);
    srf_0(3) = nrbextrude(srf_0(3), [0, 0, epsilon]);
    srf_0(3) = nrbkntins(srf_0(3), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), ...
                                    (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2), []});
    
    srf_0(4) = nrb4surf([0, 1 - epsilon], [0.5 - epsilon / 2, 1 - epsilon], [0, 1], [0.5 - epsilon / 2, 1]);
    srf_0(4) = nrbextrude(srf_0(4), [0, 0, epsilon]);
    srf_0(4) = nrbkntins(srf_0(4), {(0.5 - extension_factor * epsilon / 2) / (0.5 - epsilon / 2), 0.5, []});
    
    srf_0(5) = nrb4surf([0.5 - epsilon / 2, 1 - epsilon], [0.5 + epsilon / 2, 1 - epsilon], ...
                        [0.5 - epsilon / 2, 1], [0.5 + epsilon / 2, 1]);
    srf_0(5) = nrbextrude(srf_0(5), [0, 0, epsilon]);
    srf_0(5) = nrbkntins(srf_0(5), {0.5, 0.5, []});
    
    srf_0(6) = nrb4surf([0.5 + epsilon / 2, 1 - epsilon], [1, 1 - epsilon], [0.5 + epsilon / 2, 1], [1, 1]);
    srf_0(6) = nrbextrude(srf_0(6), [0, 0, epsilon]);
    srf_0(6) = nrbkntins(srf_0(6), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), 0.5, []});
    
    srf_0(7:12) = srf_0(1:6);
    
    for ii = 7:12
        srf_0(ii).coefs(3,:,:,:) = srf_0(ii).coefs(3,:,:,:) * (1 / epsilon - 1) + epsilon;
        srf_0(ii) = nrbkntins(srf_0(ii), {[], [], (extension_factor * epsilon / 2) / (1 - epsilon / 2)});
    end
    
    srf_0(1) = nrbkntins(srf_0(1), {[], [], 0.5});
    srf_0(2) = nrbkntins(srf_0(2), {[], [], 0.5});
    srf_0(3) = nrbkntins(srf_0(3), {[], [], 0.5});
    srf_0(4) = nrbkntins(srf_0(4), {[], [], 0.5});
    srf_0(5) = nrbkntins(srf_0(5), {[], [], 0.5});
    srf_0(6) = nrbkntins(srf_0(6), {[], [], 0.5});
    
    srf = srf_0([1:4, 6:12]);
    srf_F = srf_0(5);
end

function [problem_data, problem_data_0] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 5:7 9:18 20 22 23 25:34];
    problem_data.drchlt_sides = [2 4 8 19 21 24];

    % Simplified problem
    problem_data_0.nmnn_sides = [1 3 5 6 8:17 19 21 22 24:32]; 
    problem_data_0.drchlt_sides = [2 4 7 18 20 23]; 
end

function defeaturing_data = set_defeaturing_data(problem_data)
    defeaturing_data.omega_star_patches_in_omega0 = [1:4 6:12];

    defeaturing_data.gamma_n_sides_in_omega = [5 11 14 30];
    defeaturing_data.gamma_n_sides_in_omega0.patch = [2 4 6 11];
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = [4 2 1 5];

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data);
end
