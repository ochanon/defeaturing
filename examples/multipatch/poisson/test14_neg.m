clear
close all
clc

%% Set problem data
distances_from_dirichlet_boundary = 0.005; % 1e-2 ./ 2.^(-3:6); 
feature_radius = 0.1; 

filename = 'results/test14_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) ones(size(x)); 
problem_data.h = @(x, y, ind) zeros(size(x));

[problem_data, problem_data_0] = determineBC(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32];
method_data.nquad = [5 5];


%% Main
number_of_iterations = numel(distances_from_dirichlet_boundary);
error_h1s = zeros(1, number_of_iterations);
error_h1s_from_boundary = zeros(1, number_of_iterations);
estimator = zeros(1, number_of_iterations);
measure_of_gamma = zeros(1, number_of_iterations);

for iter = 1:number_of_iterations
    distance = distances_from_dirichlet_boundary(iter);
    fprintf('----- Distance from the Dirichlet boundary = %f -----\n', distance);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(distance, feature_radius);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    error_h1s(iter) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    [estimator(iter), measure_of_gamma(iter), error_h1s_from_boundary(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end


%% Display and save the results
if saveIt
    save(filename, 'distances_from_dirichlet_boundary', 'error_h1s', ...
        'error_h1s_from_boundary', 'estimator', 'measure_of_gamma')
end
if plotIt
    fig = figure;
    loglog(distances_from_dirichlet_boundary, error_h1s, '+-r', ...
           distances_from_dirichlet_boundary, estimator, '+-b');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', 'Location', 'west')
    xlabel('Distance between the feature and \Gamma_D')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end


%% Auxiliary functions
function [srf_0, srf, srf_F] = buildGeometry(distance, feature_radius)
    L = sqrt(2) * feature_radius / 4;
    extension_factor = 4; 

    srf_0(1) = nrbdegelev(nrbruled (nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], ...
                                    pi/4, 3*pi/4), nrbline([0.5, 0.5], [-0.5, 0.5])), [0, 1]);
    srf_0(2) = nrbdegelev(nrbruled (nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], ...
                                    3*pi/4, 5*pi/4), nrbline([-0.5, 0.5], [-0.5, -0.5])), [0, 1]);
    srf_0(3) = nrbdegelev(nrbruled (nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], ...
                                    5*pi/4, 7*pi/4), nrbline([-0.5, -0.5], [0.5, -0.5])), [0, 1]);
    srf_0(4) = nrbdegelev(nrbruled (nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], ...
                                    7*pi/4, 9*pi/4), nrbline([0.5, -0.5], [0.5, 0.5])), [0, 1]);
    for ii = 1:4
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (extension_factor - 1) * feature_radius / (1 - feature_radius)});
    end
    srf_0(5) = nrbruled (nrbline([L, -0.5 + feature_radius + distance + L], [-L, -0.5 + feature_radius + distance + L]), ...
                         nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], pi/4, 3*pi/4));
    srf_0(6) = nrbruled (nrbline([-L, -0.5 + feature_radius + distance + L], [-L, -0.5 + feature_radius + distance - L]), ...
                         nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], 3*pi/4, 5*pi/4));
    srf_0(7) = nrbruled (nrbline([-L, -0.5 + feature_radius + distance - L], [L, -0.5 + feature_radius + distance - L]), ...
                         nrbcirc(feature_radius, [0, -0.5+feature_radius+distance], 5*pi/4, 7*pi/4));
    srf_0(8) = nrbruled (nrbline([L, -0.5 + feature_radius + distance - L], [L, -0.5 + feature_radius + distance + L]), ...
                         nrbcirc(feature_radius, [0, -0.5 + feature_radius + distance], 7*pi/4, 9*pi/4));
    for ii = 5:8
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L, -0.5 + feature_radius + distance - L], [L, -0.5 + feature_radius + distance - L], ...
                                   [-L, -0.5 + feature_radius + distance + L], [L, -0.5 + feature_radius + distance + L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:4);
    srf_F = srf_0(5:9);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 5 7];
    problem_data.drchlt_sides = [2 4 6 8];
    problem_data.gamma_sides = [1 3 5 7];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
    problem_data_0.omega_patches = 1:4;
    problem_data_0.gamma_sides = cell(9, 1); 
    problem_data_0.gamma_sides(1:4) = {3, 3, 3, 3}; % relative to each patch
end
