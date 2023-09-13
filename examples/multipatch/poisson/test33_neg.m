clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.2.2, geometry \Omega_\varepsilon^2 --

%% Set problem data
epsilon_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test33_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.f = @(x, y, z) 10 * cos(3 * pi * x) .* sin(5 * pi * y) .* sin(7 * pi * z); 
problem_data.h = @(x, y, z, ind) zeros(size(x));
[problem_data, problem_data_0] = set_boundary_conditions(problem_data);

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
    [estimator(iter), measure_of_gamma(iter), error_H1s_boundary_representation(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);

    % 4b) COMPUTE THE DEFEATURING ERROR
    error_H1s(iter) = defeaturing_error_H1s(msh_0, space_0, u_0, problem_data_0.omega_patches, ...
                                            msh, space, u);

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
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', epsilon_values, epsilon_values.^(5/2), 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^{5/2}', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(epsilon)
    extension_factor = 4;

    srf_0(1) = nrb4surf ([0, 1 - epsilon], [1 - epsilon, 1 - epsilon], [0, 1], [1 - epsilon, 1]);
    srf_0(1) = nrbextrude(srf_0(1), [0, 0, epsilon]);
    srf_0(1) = nrbkntins(srf_0(1), {1 - extension_factor * epsilon / (2 * (1 - epsilon)), 0.5, []});
    
    srf_0(2) = nrb4surf ([1 - epsilon, 1 - epsilon], [1, 1 - epsilon], [1 - epsilon, 1], [1, 1]);
    srf_0(2) = nrbextrude(srf_0(2), [0, 0, epsilon]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, 0.5, []});
    
    srf_0(3) = nrb4surf ([0, 0], [1 - epsilon, 0], [0, 1 - epsilon], [1 - epsilon, 1 - epsilon]);
    srf_0(3) = nrbextrude(srf_0(3), [0, 0, epsilon]);
    srf_0(3) = nrbkntins(srf_0(3), {1 - extension_factor * epsilon / (2 * (1 - epsilon)), ...
                                    1 - extension_factor * epsilon / (2 * (1 - epsilon)), []});
    
    srf_0(4) = nrb4surf ([1 - epsilon, 0], [1, 0], [1 - epsilon, 1 - epsilon], [1, 1 - epsilon]);
    srf_0(4) = nrbextrude(srf_0(4), [0, 0, epsilon]);
    srf_0(4) = nrbkntins(srf_0(4), {0.5, 1 - extension_factor * epsilon / (2 * (1 - epsilon)), []});
    
    srf_0(5:8) = srf_0(1:4);
    
    for ii = 5:8
        srf_0(ii).coefs(3,:,:,:) = srf_0(ii).coefs(3,:,:,:) * (1 / epsilon - 1) + epsilon;
        srf_0(ii) = nrbkntins(srf_0(ii), {[], [], extension_factor * epsilon / (2 * (1 - epsilon))});
    end
    
    for ii = 1:4
        srf_0(ii) = nrbkntins(srf_0(ii), {[], [], 0.5});
    end
    
    srf = srf_0([1, 3:8]);
    srf_F = srf_0(2);
end

function [problem_data, problem_data_0] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1:5 7 8 10:19 21 22 24];
    problem_data.drchlt_sides = [6 9 20 23];
    problem_data.gamma_sides = [2 10 17]; 
    problem_data.omega0_patches = 1:7;

    % Simplified problem
    problem_data_0.nmnn_sides = [1:7 9 10 12:19 21 22 24]; 
    problem_data_0.drchlt_sides = [8 11 20 23]; 
    problem_data_0.omega_patches = [1, 3:8];
    problem_data_0.gamma_sides = cell(8,1); 
    problem_data_0.gamma_sides([1 4 6]) = {2, 4, 5}; % relative to each patch
end
