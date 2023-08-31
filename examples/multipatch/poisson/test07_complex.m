clear
close all
clc

%% Set problem data
% -- To reproduce [BCV2022], Section 6.2.1, geometry \Omega_\varepsilon^5 --
eps_values = 1e-2 ./ 2.^(0:6);
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

[problem_data, problem_data_0, problem_data_F] = determineBC(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32]; 
method_data.nquad = [5 5];


%% Main
number_of_epsilons = numel(eps_values);
error_h1s = zeros(1, number_of_epsilons);
error_h1s_Omega_star = zeros(1, number_of_epsilons);
error_h1s_F = zeros(1, number_of_epsilons);
error_h1s_from_boundary = zeros(1, number_of_epsilons);
error_h1s_from_positive_boundary = zeros(1, number_of_epsilons);
error_h1s_from_negative_boundary = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
estimator_positive_feature = zeros(1, number_of_epsilons);
estimator_negative_feature = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_h1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = eps_values(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F_positive] = buildGeometry(epsilon); 
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_F.geo_name = srf_F_positive;

    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);

    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    problem_data_F.h = buildDirichletConditionExtPb(msh_0, problem_data_0, space_0, u_0);

    % 3b) SOLVE THE EXTENSION PROBLEM
    [F, msh_F, space_F, u_0tilde] = my_mp_solve_laplace (problem_data_F, method_data);

    % 4a) COMPUTE ERROR
    [error_h1s(iter), error_h1s_Omega_star(iter), error_h1s_F(iter)] = errh1s_positive(msh, space, u, msh_0, ...
        space_0, u_0, msh_F, space_F, u_0tilde, problem_data.omega0_patches, ...
        problem_data_F.F_patches, problem_data_0.omega_patches);

    % 4b) COMPUTE ESTIMATOR AND ERROR FROM BOUNDARY
    [estimator_positive_feature(iter), ~, ~, error_h1s_from_positive_boundary(iter)] = ...
        est_positive(msh_F, space_F, u_0tilde, problem_data_F.gamma0_sides, problem_data_F.gammae_sides,...
            problem_data_0.g, problem_data.g, problem_data_F.F_patches, ...
            problem_data.omega0_patches, space, u);
    [estimator_negative_feature(iter), ~, error_h1s_from_negative_boundary(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
    
    norm_of_u(iter) = errh1s_positive(msh, space, u, ...
        msh_0, space_0, zeros(size(u_0)), msh_F, space_F, zeros(size(u_0tilde)), problem_data.omega0_patches, ...
        problem_data_F.F_patches, problem_data_0.omega_patches);
    relative_error_h1s(iter) = error_h1s(iter)/norm_of_u(iter);

    estimator(iter) = sqrt(estimator_positive_feature(iter)^2 + estimator_negative_feature(iter)^2); 
    error_h1s_from_boundary(iter) = sqrt(error_h1s_from_positive_boundary(iter)^2 ...
                                          + error_h1s_from_negative_boundary(iter)^2);
end


%% Display and save the results
if saveIt
    save(filename, 'eps_values', 'error_h1s', 'error_h1s_Omega_star', 'error_h1s_F', ...
        'error_h1s_from_boundary', 'error_h1s_from_positive_boundary', 'error_h1s_from_negative_boundary',...
        'estimator', 'estimator_positive_feature', 'estimator_negative_feature', 'norm_of_u', 'relative_error_h1s')
end
if plotIt
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, estimator, '+-b', ...
           eps_values, eps_values, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf("For epsilon (feature's side length) = %e, \n", epsilon)
fprintf('    * Estimator                           = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_h1s(end))
fprintf('    * Effectivity index                   = %f \n', estimator(end)/error_h1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_F_positive] = buildGeometry(epsilon)
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

function [problem_data, problem_data_0, problem_data_F] = determineBC(problem_data)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F = rmfield(problem_data_F, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 5 6 8:16]; 
    problem_data.drchlt_sides = [2 3 4 7];
    problem_data.gamma_sides = [5 10 11]; % for the negative feature
    problem_data.omega0_patches = 1:7;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 5 7:12];
    problem_data_0.drchlt_sides = [2 3 4 6];
    problem_data_0.gamma0_sides = 9; % for the Dirichlet extension
    problem_data_0.omega_patches = [1:6 8];
    problem_data_0.gamma_sides = cell(8, 1); % relative to each patch
    problem_data_0.gamma_sides([3 6 8]) = {4, 2, 1};
    
    % Extension problem
    problem_data_F.nmnn_sides = [1 2 4]; 
    problem_data_F.drchlt_sides = 3; 
    problem_data_F.gamma0_sides = 3; 
    problem_data_F.gammae_sides = {[]}; % relative to each patch
    problem_data_F.F_patches = 1;
end
