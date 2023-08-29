clear
close all
clc

%% Set problem data
eps_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test04_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 40 * cos(pi * x) .* cos(pi * y) + 10 * cos(5 * pi * x) .* cos(7 * pi * y);
[problem_data, problem_data_0] = determineBC(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32]; 
method_data.nquad = [5 5];


%% Main
number_of_epsilons = numel(eps_values);
error_h1s = zeros(1, number_of_epsilons);
error_h1s_from_interface = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_h1s = zeros(1, number_of_epsilons);
measure_of_gamma = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = eps_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    error_h1s(iter) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    norm_of_u(iter) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    relative_error_h1s(iter) = error_h1s(iter)/norm_of_u(iter);
 
    [estimator(iter), measure_of_gamma(iter), error_h1s_from_interface(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end


%% Display and save the results
if saveIt
    save(filename, 'eps_values', 'error_h1s', 'error_h1s_from_interface', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_h1s')
end
if plotIt
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, error_h1s_from_interface, '+-.r', ...
           eps_values, estimator, '+-b', eps_values, eps_values*18, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', '\int_\gamma d_\gamma (u-u_0)', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = buildGeometry(epsilon)
    extension_factor = 4;

    srf_0(1) = nrb4surf([0, 0], [0.5 - epsilon / 2, 0], [0, 1 - epsilon], [0.5 - epsilon / 2, 1 - epsilon]);
    srf_0(1) = nrbkntins(srf_0(1), {(0.5 - extension_factor * epsilon / 2) / (0.5 - epsilon / 2), ...
                                    (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2)});

    srf_0(2) = nrb4surf([0.5 - epsilon/2, 0], [0.5 + epsilon / 2, 0], ...
                        [0.5 - epsilon / 2, 1 - epsilon], [0.5 + epsilon / 2, 1 - epsilon]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2)});

    srf_0(3) = nrb4surf([0.5 + epsilon / 2, 0], [1, 0], [0.5 + epsilon / 2, 1 - epsilon], [1, 1-epsilon]);
    srf_0(3) = nrbkntins(srf_0(3), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), ...
                                    (1 - extension_factor * epsilon / 2) / (1 - epsilon / 2)});

    srf_0(4) = nrb4surf([0, 1 - epsilon], [0.5 - epsilon / 2, 1 - epsilon], [0, 1], [0.5 - epsilon / 2, 1]);
    srf_0(4) = nrbkntins(srf_0(4), {(0.5 - extension_factor * epsilon / 2) / (0.5 - epsilon / 2), 0.5});

    srf_0(5) = nrb4surf([0.5 - epsilon / 2, 1 - epsilon], [0.5 + epsilon / 2, 1 - epsilon], ...
                        [0.5 - epsilon / 2, 1], [0.5 + epsilon / 2, 1]);
    srf_0(5) = nrbkntins(srf_0(5), {0.5, 0.5});

    srf_0(6) = nrb4surf([0.5 + epsilon / 2, 1 - epsilon], [1, 1 - epsilon], [0.5 + epsilon / 2, 1], [1, 1]);
    srf_0(6) = nrbkntins(srf_0(6), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), 0.5});

    srf = srf_0([1:4, 6]);
    srf_F = srf_0(5);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 4 5 7 8 9 10 11 12];
    problem_data.drchlt_sides = [2 3 6];
    problem_data.gamma_sides = [4 8 10];
    problem_data.omega0_patches = 1:5;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 4 6 7 8 9 10];
    problem_data_0.drchlt_sides = [2 3 5];
    problem_data_0.omega_patches = [1:4 6];
    problem_data_0.gamma_sides = cell(6, 1); 
    problem_data_0.gamma_sides([2 4 6]) = {4, 2, 1}; % relative to each patch
end