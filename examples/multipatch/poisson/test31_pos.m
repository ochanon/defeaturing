clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.2.2, geometry \Omega_\varepsilon^3 --

%% Set problem data
epsilon_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test31_pos';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y, z) ones(size(x));
problem_data.g = @(x, y, z, ind) zeros(size(x));
problem_data.f = @(x, y, z) 10 * cos(3 * pi * x) .* sin(5 * pi * y) .* sin(7 * pi * z); 
problem_data.h = @(x, y, z, ind) zeros(size(x));

[problem_data, problem_data_0, problem_data_F] = set_boundary_conditions(problem_data);
defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_F);

method_data.degree = [3 3 3];
method_data.regularity = [2 2 2];
method_data.nsub = [8 8 8];
method_data.nquad = [5 5 5];


%% Main
number_of_epsilons = numel(epsilon_values);
error_H1s = zeros(1, number_of_epsilons);
error_H1s_0 = zeros(1, number_of_epsilons);
error_H1s_F = zeros(1, number_of_epsilons);
error_H1s_boundary_representation = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
measure_of_gamma_0 = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_H1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = epsilon_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = build_geometry(epsilon);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_F.geo_name = srf_F;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_generalized(problem_data, method_data);
    
    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_generalized(problem_data_0, method_data);
    problem_data_F.h = extract_boundary_dofs(msh_0, space_0, u_0, defeaturing_data.gamma_0p_sides_in_omega0);
    
    % 3b) SOLVE THE EXTENSION PROBLEM
    [F, msh_F, space_F, u_0tilde] = mp_solve_laplace_generalized(problem_data_F, method_data);
    
    % 4a) COMPUTE THE DEFEATURING ESTIMATOR
    [estimator(iter), measure_of_gamma_0(iter)] = ...
        estimate_defeaturing_error_H1s(defeaturing_data, msh_0, space_0, u_0, msh_F, space_F, u_0tilde);

    % 4b) COMPUTE THE DEFEATURING ERROR
    [error_H1s(iter), error_H1s_0(iter), error_H1s_F(iter)] = defeaturing_error_H1s(msh_0, space_0, u_0, [], ...
                                                  msh, space, u, defeaturing_data.omega_star_patches_in_omega, ...
                                                  msh_F, space_F, u_0tilde);
    error_H1s_boundary_representation(iter) = ...
        defeaturing_error_H1s_boundary_representation(defeaturing_data, msh, space, u, ...
                                                      msh_0, space_0, u_0, msh_F, space_F, u_0tilde);

    norm_of_u(iter) = error_H1s_in_patches(msh, space, u, 1:msh.npatch, msh, space, zeros(size(u)));
    relative_error_H1s(iter) = error_H1s(iter) / norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'epsilon_values', 'error_H1s', 'error_H1s_0', 'error_H1s_F', ...
        'error_H1s_boundary_representation', 'estimator', 'measure_of_gamma_0', ...
        'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', ...
           epsilon_values, epsilon_values.^(3/2)*1e-1, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^{3/2}', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
    
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, error_H1s_0, '+-g', ...
           epsilon_values, error_H1s_F, '+-c', epsilon_values, estimator, '+-b', ...
           epsilon_values, epsilon_values.^(3/2)*1e-1, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', '|u-u_0|_{1,\Omega_0}', '|u-u_0|_{1,F}',...
        'Estimator', '\epsilon^{3/2}', 'Location', 'northwest')
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(epsilon)
    extension_factor = 4;

    srf(1) = nrb4surf ([0, 0], [0.5 - epsilon / 2, 0], [0, 1], [0.5 - epsilon / 2, 1]);
    srf(1) = nrbextrude(srf(1), [0, 0, epsilon]);
    srf(1) = nrbkntins(srf(1), {0.5, 1 - extension_factor * epsilon / 2, []}); 
    
    srf(2) = nrb4surf ([0.5 - epsilon / 2, 0], [0.5 + epsilon / 2, 0], [0.5 - epsilon / 2, 1], [0.5 + epsilon / 2, 1]);
    srf(2) = nrbextrude(srf(2), [0, 0, epsilon]);
    srf(2) = nrbkntins(srf(2), {0.5, 1 - extension_factor * epsilon / 2, []});
    
    srf(3) = nrb4surf ([0.5 + epsilon / 2, 0], [1, 0], [0.5 + epsilon / 2, 1], [1, 1]);
    srf(3) = nrbextrude(srf(3), [0, 0, epsilon]);
    srf(3) = nrbkntins(srf(3), {0.5, 1 - extension_factor * epsilon / 2, []}); 
    
    srf(4:6) = srf(1:3);
    
    for ii = 4:6
        srf(ii).coefs(3,:,:,:) = srf(ii).coefs(3,:,:,:) * (1 / epsilon - 1) + epsilon;
        srf(ii) = nrbkntins(srf(ii), {[], [], ((extension_factor - 1) * epsilon / 2) / (1 - epsilon / 2)});
    end
    
    srf(7) = nrb4surf([0.5 - epsilon / 2, 1], [0.5 + epsilon / 2, 1], ...
                      [0.5 - epsilon / 2, 1 + epsilon], [0.5 + epsilon / 2, 1 + epsilon]);
    srf(7) = nrbextrude(srf(7), [0, 0, epsilon]);
    srf(7) = nrbkntins(srf(7), {0.5, 0.5, 0.5});
    
    srf(1) = nrbkntins(srf(1), {[], [], 0.5});
    srf(2) = nrbkntins(srf(2), {[], [], 0.5});
    srf(3) = nrbkntins(srf(3), {[], [], 0.5});
    
    srf_0 = srf(1:6);
    srf_F = srf(7);
end

function [problem_data, problem_data_0, problem_data_F] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F = rmfield(problem_data_F, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 4 6 7 9:11 13 14 16:18 20:26]; 
    problem_data.drchlt_sides = [2 5 8 12 15 19];

    % Simplified problem
    problem_data_0.nmnn_sides = [1 3 4 6:8 10:12 14 15 17:19 21 22];
    problem_data_0.drchlt_sides = [2 5 9 13 16 20];
    
    % Extension problem
    problem_data_F.nmnn_sides = [1 2 4:6]; 
    problem_data_F.drchlt_sides = 3; 
end

function defeaturing_data = set_defeaturing_data(problem_data, problem_data_0, problem_data_F)
    defeaturing_data.omega_star_patches_in_omega = 1:6;
    defeaturing_data.Fp_patches_in_tildeFp = 1;

    defeaturing_data.gamma_0p_sides_in_omega0 = 6; 
    defeaturing_data.gamma_0p_sides_in_tildeFp = 3; 

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data, problem_data_0, problem_data_F);
end
