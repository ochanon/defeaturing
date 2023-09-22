clear
close all
clc

%% Set problem data
% -- To reproduce [BCV2022], Section 6.1.1, geometry \Omega_c --
%  * Case r_c = 6.37e−2, uncomment the three lines below:
% perimeter = 0.4;
% circle_radiuses = perimeter/(2*pi);
% neumann_function = @(x, y, ind, radius) zeros(size(x));
% --
%  * Case r_c = 5.64e−2, uncomment the three lines below:
% area = 1e-2;
% circle_radiuses = sqrt(area/pi);
% neumann_function = @(x, y, ind, radius) zeros(size(x));
% --

% -- To reproduce [BCV2022], Section 6.2.3, geometry \Omega_\varepsilon --
circle_radiuses = 1e-2 ./ 2.^(0:6);
g1 = @(x, y, ind, radius) zeros(size(x));
g2 = @(x, y, ind, radius) ones(size(x));
g3 = @(x, y, ind, radius) ones(size(x)) ./ radius;
g4 = @(x, y, ind, radius) ones(size(x)) ./ (radius.^3);
neumann_function = g4; % Choose here the required Neumann data g1, g2, g3 or g4
% --

filename = 'results/test18_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.f = @(x, y) ones(size(x)); 
problem_data.h = @(x, y, ind) zeros(size(x));

[problem_data, problem_data_0] = set_boundary_conditions(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32];
method_data.nquad = [5 5];


%% Main
number_of_circle_radiuses = numel(circle_radiuses);
error_H1s = zeros(1, number_of_circle_radiuses);
error_H1s_boundary_representation = zeros(1, number_of_circle_radiuses);
estimator = zeros(1, number_of_circle_radiuses);
measure_of_gamma = zeros(1, number_of_circle_radiuses);
norm_of_u = zeros(1, number_of_circle_radiuses);
relative_error_H1s = zeros(1, number_of_circle_radiuses);

for iter = 1:number_of_circle_radiuses
    radius = circle_radiuses(iter);
    fprintf('----- circle radius = %f -----\n', radius);
    problem_data.g = @(x, y, ind) neumann_function(x, y, ind, radius);
    defeaturing_data = set_defeaturing_data(problem_data);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = build_geometry(radius);
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
    save(filename, 'circle_radiuses', 'error_H1s', 'error_H1s_boundary_representation', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(circle_radiuses, error_H1s, '+-r', circle_radiuses, estimator, '+-b');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf('For circle_radius = %e, \n', circle_radiuses(end))
fprintf('    * Estimator                                      = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}            = %e \n', error_H1s(end))
fprintf('    * Relative error |u-u_0|_{1,Omega}/|u|_{1,Omega} = %e \n', relative_error_H1s(end))
fprintf('    * Effectivity index                              = %f \n', estimator(end) / error_H1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(epsilon)
    L = sqrt(2) * epsilon / 4;
    extension_factor = 4;

    srf_0(1) = nrbruled (nrbcirc(epsilon, [0 0], pi/4, 3*pi/4), nrbcirc(1, [0 0], pi/4, 3*pi/4));
    srf_0(2) = nrbtform(srf_0(1), vecrotz(pi/2));
    srf_0(3) = nrbtform(srf_0(2), vecrotz(pi/2));
    srf_0(4) = nrbtform(srf_0(3), vecrotz(pi/2));
    for ii = 1:4
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (extension_factor-1) * epsilon / (1-epsilon)});
    end
    srf_0(5) = nrbruled (nrbline([L L], [-L L]), nrbcirc(epsilon, [0 0], pi/4, 3*pi/4));
    srf_0(6) = nrbtform(srf_0(5), vecrotz(pi/2));
    srf_0(7) = nrbtform(srf_0(6), vecrotz(pi/2));
    srf_0(8) = nrbtform(srf_0(7), vecrotz(pi/2));
    for ii = 5:8
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L -L], [L, -L], [-L L], [L L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:4);
    srf_F = srf_0(5:9);
end

function [problem_data, problem_data_0] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 5 7];
    problem_data.drchlt_sides = [2 4 6 8];

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
end

function defeaturing_data = set_defeaturing_data(problem_data)
    defeaturing_data.omega_star_patches_in_omega0 = 1:4;

    defeaturing_data.gamma_n_sides_in_omega = [1 3 5 7];
    defeaturing_data.gamma_n_sides_in_omega0.patch = 1:4;
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = [3 3 3 3];

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data);
end
