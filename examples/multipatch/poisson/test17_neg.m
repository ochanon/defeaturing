clear
close all
clc

%% Set problem data
% -- To reproduce [BCV2022], Section 6.1.1, geometry \Omega_s --
%  * Case r_s = 5e-2:
half_side_length_of_squared_feature = 5 * 1e-2;
% --

% -- To study the convergence with respect to the size of the feature --
% half_side_length_of_squared_feature = 1e-2 ./ 2.^(1:7);
% --

side_lengths = 2 * half_side_length_of_squared_feature;
filename = 'results/test17_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) ones(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

[problem_data, problem_data_0] = set_boundary_conditions(problem_data);
defeaturing_data = set_defeaturing_data(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [16 16];
method_data.nquad = [5 5];


%% Main
number_of_side_lengths = numel(side_lengths);
error_H1s = zeros(1, number_of_side_lengths);
error_H1s_boundary_representation = zeros(1, number_of_side_lengths);
estimator = zeros(1, number_of_side_lengths);
measure_of_gamma = zeros(1, number_of_side_lengths);
norm_of_u = zeros(1, number_of_side_lengths);
relative_error_H1s = zeros(1, number_of_side_lengths);

for iter = 1:number_of_side_lengths
    side = side_lengths(iter);
    fprintf('----- side length = %f -----\n', side);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = build_geometry(side);
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
    save(filename, 'side_lengths', 'error_H1s', 'error_H1s_boundary_representation', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(side_lengths, error_H1s, '+-r', side_lengths, estimator, '+-b', ...
           side_lengths, side_lengths.^2*2, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^2', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf('For half_side_length_of_squared_feature = %e, \n', half_side_length_of_squared_feature(end))
fprintf('    * Estimator                                      = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}            = %e \n', error_H1s(end))
fprintf('    * Relative error |u-u_0|_{1,Omega}/|u|_{1,Omega} = %e \n', relative_error_H1s(end))
fprintf('    * Effectivity index                              = %f \n', estimator(end) / error_H1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(side_length)
    L = side_length/2;
    extension_factor = 4;
    tol = 0.1;

    srf_0(1) = nrbruled (nrbcirc(tol, [0 0], pi/4, 3*pi/4), nrbcirc(1, [0 0], pi/4, 3*pi/4));
    srf_0(2) = nrbtform(srf_0(1), vecrotz(pi/2));
    srf_0(3) = nrbtform(srf_0(2), vecrotz(pi/2));
    srf_0(4) = nrbtform(srf_0(3), vecrotz(pi/2));
    
    srf_0(5) = nrbruled (nrbline([L L], [-L L]), nrbcirc(tol, [0 0], pi/4, 3*pi/4));
    srf_0(6) = nrbtform(srf_0(5), vecrotz(pi/2));
    srf_0(7) = nrbtform(srf_0(6), vecrotz(pi/2));
    srf_0(8) = nrbtform(srf_0(7), vecrotz(pi/2));

    if extension_factor * side_length < 0.1
        for ii = 1:4
            srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
            srf_0(ii+4) = nrbkntins(srf_0(ii+4), {0.5, extension_factor * side_length / 0.1});
        end
    else
        for ii = 1:4
            srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (extension_factor * side_length - 0.1) / 0.9});
            srf_0(ii+4) = nrbkntins(srf_0(ii + 4), {0.5, 0.5});
        end
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L -L], [L, -L], [-L L], [L L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:8);
    srf_F = srf_0(9);
end

function [problem_data, problem_data_0] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = 5:8;
    problem_data.drchlt_sides = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
end

function defeaturing_data = set_defeaturing_data(problem_data)
    defeaturing_data.omega_star_patches_in_omega0 = 1:8;

    defeaturing_data.gamma_n_sides_in_omega = 5:8;
    defeaturing_data.gamma_n_sides_in_omega0.patch = 5:8;
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = [3 3 3 3];

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data);
end
