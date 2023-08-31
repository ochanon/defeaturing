clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.2.1, geometry \Omega_\varepsilon^1 --

%% Set problem data
eps_values = 1e-2 ./ 2.^(0:6);

filename = 'results/test16_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 10 * cos(3 * pi * x) .* sin(5 * pi * y); 
problem_data.h = @(x, y, ind) zeros(size(x));
[problem_data, problem_data_0] = determineBC(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32];
method_data.nquad = [5 5];


%% Main
number_of_epsilons = numel(eps_values);
error_h1s = zeros(1, number_of_epsilons);
error_h1s_from_boundary = zeros(1, number_of_epsilons);
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
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    error_h1s(iter) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    norm_of_u(iter) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    relative_error_h1s(iter) = error_h1s(iter)/norm_of_u(iter);

    [estimator(iter), measure_of_gamma(iter), error_h1s_from_boundary(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end


%% Display and save the results
if saveIt
    save(filename, 'eps_values', 'error_h1s', 'error_h1s_from_boundary', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_h1s')
end
if plotIt
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, estimator, '+-b', eps_values, eps_values, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = buildGeometry(epsilon)
    s = epsilon/(2*sqrt(2));
    extension_factor = 4;

    srf_0(1) = nrbdegelev(nrbsquare([0.5-s, 1-s], 2*s, s), [1,1]);
    srf_0(1) = nrbkntins(srf_0(1), {0.5, 0.5});
    
    crv1 = nrbline([0.5-s, 1], [0.5-s, 1-s]);
    crv2 = nrbcirc(epsilon, [0.5, 1], pi, 5*pi/4);
    srf_0(2) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, 0.5});
    
    crv1 = nrbline([0.5-s, 1-s], [0.5+s, 1-s]);
    crv2 = nrbcirc(epsilon, [0.5, 1], 5*pi/4, 7*pi/4);
    srf_0(3) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(3) = nrbkntins(srf_0(3), {0.5, 0.5});
    
    crv1 = nrbline([0.5+s, 1-s], [0.5+s, 1]);
    crv2 = nrbcirc(epsilon, [0.5, 1], -pi/4, 0);
    srf_0(4) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(4) = nrbkntins(srf_0(4), {0.5, 0.5});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], pi, 5*pi/4);
    crv2 = nrbline([0, 1], [0, 0]);
    srf_0(5) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(5) = nrbkntins(srf_0(5), {0.5, (extension_factor - 1) * epsilon / (1 - epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], 5*pi/4, 7*pi/4);
    crv2 = nrbline([0, 0], [1, 0]);
    srf_0(6) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(6) = nrbkntins(srf_0(6), {0.5, (extension_factor - 1) * epsilon / (1 - epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], -pi/4, 0);
    crv2 = nrbline([1, 0], [1, 1]);
    srf_0(7) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(7) = nrbkntins(srf_0(7), {0.5, (extension_factor - 1) * epsilon / (1 - epsilon)});
    
    srf_F = srf_0(1:4);
    srf = srf_0(5:7);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 2 4 6 7];
    problem_data.drchlt_sides = [3 5 8];
    problem_data.gamma_sides = [2 4 7];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 2 3 4 7];
    problem_data_0.drchlt_sides = [5 6 8];
    problem_data_0.omega_patches = 5:7;
    problem_data_0.gamma_sides = cell(7,1); 
    problem_data_0.gamma_sides(5:7) = {3, 3, 3}; % relative to each patch
end
