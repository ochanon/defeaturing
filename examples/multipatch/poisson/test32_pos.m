clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.2.2, geometry \Omega_\varepsilon^4 --

%% Set problem data
eps_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test32_pos';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.f = @(x, y, z) 10 * cos(3 * pi * x) .* sin(5 * pi * y) .* sin(7 * pi * z); 
problem_data.h = @(x, y, z, ind) zeros(size(x));
[problem_data, problem_data_0, problem_data_ext] = determineBC(problem_data);

method_data.degree = [3 3 3];
method_data.regularity = [2 2 2];
method_data.nsub = [8 8 8];
method_data.nquad = [5 5 5];


%% Main
number_of_epsilons = numel(eps_values);
error_h1s = zeros(1, number_of_epsilons);
error_h1s_0 = zeros(1, number_of_epsilons);
error_h1s_F = zeros(1, number_of_epsilons);
error_h1s_from_interface = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
measure_of_gamma = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_h1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = eps_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_ext.geo_name = srf_F;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    problem_data_ext.h = buildDirichletConditionExtPb(msh_0, problem_data_0, space_0, u_0);
    % 3b) SOLVE THE EXTENSION PROBLEM
    [F, msh_F, space_F, u_0tilde] = my_mp_solve_laplace (problem_data_ext, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    [error_h1s(iter), error_h1s_0(iter), error_h1s_F(iter)] = errh1s_positive(msh, space, u, ...
        msh_0, space_0, u_0, msh_F, space_F, u_0tilde, problem_data.omega0_patches);
    [estimator(iter), measure_of_gamma(iter), ~, error_h1s_from_interface(iter)] = ...
        est_positive(msh_F, space_F, u_0tilde, problem_data_ext.gamma0_sides, problem_data_ext.gammae_sides, ...
        problem_data_0.g, problem_data.g, problem_data_ext.F_patches, problem_data.omega0_patches, space, u);
    
    norm_of_u(iter) = errh1s_positive(msh, space, u, ...
        msh_0, space_0, zeros(size(u_0)), msh_F, space_F, zeros(size(u_0tilde)), problem_data.omega0_patches);
    relative_error_h1s(iter) = error_h1s(iter)/norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'eps_values', 'error_h1s', 'error_h1s_0', 'error_h1s_F', ...
        'error_h1s_from_interface', 'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_h1s')
end
if plotIt
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, estimator, '+-b', eps_values, eps_values.^(5/2)*0.7, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^{5/2}', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
    
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, error_h1s_0, '+-g', ...
        eps_values, error_h1s_F, '+-c', eps_values, estimator, '+-b', eps_values, eps_values.^(5/2)*0.7, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', '|u-u_0|_{1,\Omega_0}', '|u-u_0|_{1,F}',...
        'Estimator', '\epsilon^{5/2}', 'Location', 'northwest')
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = buildGeometry(epsilon)
    extension_factor = 4;

    srf(1) = nrb4surf([0, 0], [1 - epsilon, 0], [0, 1], [1 - epsilon, 1]);
    srf(1) = nrbextrude(srf(1), [0, 0, epsilon]);
    srf(1) = nrbkntins(srf(1), {(1 - extension_factor * epsilon) / (1 - epsilon), 1 - extension_factor * epsilon, []});
    
    srf(2) = nrb4surf([1 - epsilon, 0], [1, 0], [1 - epsilon, 1], [1, 1]);
    srf(2) = nrbextrude(srf(2), [0, 0, epsilon]);
    srf(2) = nrbkntins(srf(2), {0.5, 1 - extension_factor * epsilon, []});
    
    srf(3:4) = srf(1:2);
    
    for ii = 3:4
        srf(ii).coefs(3,:,:,:) = srf(ii).coefs(3,:,:,:) * (1 / epsilon - 1) + epsilon;
        srf(ii) = nrbkntins(srf(ii), {[], [], 1 - (1 - extension_factor * epsilon) / (1 - epsilon)});
    end
    
    srf(5) = nrb4surf([1 - epsilon, 1], [1, 1], [1 - epsilon, 1 + epsilon], [1, 1 + epsilon]);
    srf(5) = nrbextrude(srf(5), [0, 0, epsilon]);
    srf(5) = nrbkntins(srf(5), {0.5, 0.5, 0.5});
    
    srf(1) = nrbkntins(srf(1), {[], [], 0.5});
    srf(2) = nrbkntins(srf(2), {[], [], 0.5});
    
    srf_0 = srf(1:4);
    srf_F = srf(5);
end

function [problem_data, problem_data_0, problem_data_ext] = determineBC(problem_data)
    problem_data_0 = problem_data;
    problem_data_ext = problem_data;
    problem_data_ext = rmfield(problem_data_ext, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 3:5 7:8 10:12 14:20]; 
    problem_data.drchlt_sides = [2 6 9 13];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 3:5 7:9 11:13 15:16];
    problem_data_0.drchlt_sides = [2 6 10 14];
    problem_data_0.gamma0_sides = 7; 
    
    % Extension problem
    problem_data_ext.nmnn_sides = [1 2 4:6]; 
    problem_data_ext.drchlt_sides = 3; 
    problem_data_ext.gamma0_sides = 3; 
    problem_data_ext.gammae_sides = {[]};
    problem_data_ext.F_patches = 1;
end
