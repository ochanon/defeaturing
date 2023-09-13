clear
close all
clc

%% Set problem data
epsilon_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test09_pos';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 10 * cos(3 * pi * x) .* sin(5 * pi * y); 
problem_data.h = @(x, y, ind) zeros(size(x));
[problem_data, problem_data_0, problem_data_Ftilde] = set_boundary_conditions(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [16 16];
method_data.nquad = [5 5];


%% Main
number_of_epsilons = numel(epsilon_values);
error_H1s = zeros(1, number_of_epsilons);
error_H1s_0 = zeros(1, number_of_epsilons);
error_H1s_F = zeros(1, number_of_epsilons);
error_H1s_boundary_representation = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
measure_of_gammae = zeros(1, number_of_epsilons);
measure_of_gamma0 = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = epsilon_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F, srf_Ftilde] = build_geometry(epsilon);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_Ftilde.geo_name = srf_Ftilde;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_generalized(problem_data, method_data);
    
    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_generalized(problem_data_0, method_data);
    problem_data_Ftilde.h = extract_boundary_dofs(msh_0, space_0, u_0, problem_data_0.gamma0_sides);

    % 3b) SOLVE THE EXTENSION PROBLEM
    [Ftilde, msh_Ftilde, space_Ftilde, u_0tilde] = mp_solve_laplace_generalized(problem_data_Ftilde, method_data);
    
    % 4a) COMPUTE THE DEFEATURING ESTIMATOR
    [estimator(iter), measure_of_gamma0(iter), measure_of_gammae(iter), error_H1s_boundary_representation(iter)] = ...
        est_positive(msh_Ftilde, space_Ftilde, u_0tilde, problem_data_Ftilde.gamma0_sides,...
            problem_data_Ftilde.gammae_sides, problem_data_0.g, problem_data.g, ...
            problem_data_Ftilde.F_patches, problem_data.omega0_patches, space, u);
    
    % 4b) COMPUTE THE DEFEATURING ERROR
    [error_H1s(iter), error_H1s_0(iter), error_H1s_F(iter)] = ...
        defeaturing_error_H1s(msh_0, space_0, u_0, [], ...
                              msh, space, u, problem_data.omega0_patches, ...
                              msh_Ftilde, space_Ftilde, u_0tilde, problem_data_Ftilde.F_patches);
end


%% Display and save the results
if saveIt
    save(filename, 'epsilon_values', 'error_H1s', 'error_H1s_0', 'error_H1s_F', ...
        'error_H1s_boundary_representation', 'estimator', 'measure_of_gamma0', 'measure_of_gammae')
end
if plotIt
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', epsilon_values, epsilon_values*0.3, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
    
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, error_H1s_0, '+-g', ...
        epsilon_values, error_H1s_F, '+-c', epsilon_values, estimator, '+-b', epsilon_values, epsilon_values*0.3, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', '|u-u_0|_{1,\Omega_0}', '|u-u_0|_{1,F}',...
        'Estimator', '\epsilon', 'Location', 'northwest')
end


%% Auxiliary functions
function [srf_0, srf, srf_F, srf_Ftilde] = build_geometry(epsilon)
    extension_factor = 4; 

    srf(1) = nrb4surf([0, 0], [0.5 - epsilon / 2, 0], [0, 1], [0.5 - epsilon / 2, 1]);
    srf(1) = nrbkntins(srf(1), {(1 - extension_factor * epsilon) / (1 - epsilon), ...
                                 1 - (extension_factor - 1) * epsilon / 2});
    
    srf(2) = nrb4surf([0.5 - epsilon / 2, 0], [0.5 - epsilon / 4, 0], ...
                      [0.5 - epsilon / 2, 1], [0.5 - epsilon / 4, 1]);
    srf(2) = nrbkntins(srf(2), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf(3) = nrb4surf([0.5 - epsilon / 4, 0], [0.5 + epsilon / 2, 0], ...
                      [0.5 - epsilon / 4, 1], [0.5 + epsilon / 2, 1]);
    srf(3) = nrbkntins(srf(3), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf(4) = nrb4surf([0.5 + epsilon / 2, 0], [1, 0], [0.5 + epsilon / 2, 1], [1, 1]);
    srf(4) = nrbkntins(srf(4), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), ...
                                1 - (extension_factor - 1) * epsilon / 2});
    
    srf(5) = nrb4surf([0.5 - epsilon / 2, 1], [0.5 - epsilon / 4, 1], ...
                      [0.5 - epsilon / 2, 1 + epsilon / 4], [0.5 - epsilon / 4, 1 + epsilon / 4]);
    srf(5) = nrbkntins(srf(5), {0.5, 0.5});
    
    srf(6) = nrb4surf([0.5 - epsilon / 2, 1 + epsilon / 4], [0.5 - epsilon / 4, 1 + epsilon / 4], ...
                      [0.5 - epsilon / 2, 1 + epsilon], [0.5, 1 + epsilon / 2]);
    srf(6) = nrbkntins(srf(6), {0.5, 0.5});
    
    srf(7) = nrb4surf([0.5 - epsilon / 4, 1 + epsilon / 4], [0.5 - epsilon / 4, 1], ...
                      [0.5, 1 + epsilon / 2], [0.5 + epsilon / 2, 1]);
    srf(7) = nrbkntins(srf(7), {0.5, 0.5});
    
    srf(8) = nrb4surf([0.5 + epsilon / 4, 1 + 3 * epsilon / 4], [0.5 + epsilon / 2, 1 + 3 * epsilon / 4], ...
                      [0.5 + epsilon / 4, 1 + epsilon], [0.5 + epsilon / 2, 1 + epsilon]);
    srf(8) = nrbkntins(srf(8), {0.5, 0.5});
    
    srf(9) = nrb4surf([0.5 - epsilon / 2, 1 + epsilon], [0.5, 1 + epsilon / 2], ...
                      [0.5 + epsilon / 4, 1 + epsilon], [0.5 + epsilon / 4, 1 + 3 * epsilon / 4]);
    srf(9) = nrbkntins(srf(9), {0.5, 0.5});
    
    srf(10) = nrb4surf([0.5, 1 + epsilon / 2], [0.5 + epsilon / 4, 1 + 3 * epsilon / 4], ...
                       [0.5 + epsilon / 2, 1], [0.5 + epsilon / 2, 1 + 3 * epsilon / 4]);
    srf(10) = nrbkntins(srf(10), {0.5, 0.5});

    srf_0 = srf(1:4);
    srf_F = srf(5:7);
    srf_Ftilde = srf(5:10);
    srf = srf(1:7);
end

function [problem_data, problem_data_0, problem_data_Ftilde] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    problem_data_Ftilde = problem_data;
    problem_data_Ftilde = rmfield(problem_data_Ftilde, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 6 8:12]; 
    problem_data.drchlt_sides = [2 4 5 7];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 3 5 7 8 10];
    problem_data_0.drchlt_sides = [2 4 6 9];
    problem_data_0.gamma0_sides = [5 7]; 
    
    % Extension problem
    problem_data_Ftilde.nmnn_sides = [1 3 5:8]; 
    problem_data_Ftilde.drchlt_sides = [2 4]; 
    problem_data_Ftilde.gamma0_sides = [2 4]; 
    problem_data_Ftilde.F_patches = 1:3;
    problem_data_Ftilde.gammae_sides = {[], 4, 4};
end
