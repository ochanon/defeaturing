clear
close all
clc

%% Set problem data
number_of_star_branches = 10;
% Note: interior_radius = perimeter / (2 * number_of_star_branches 
%                                      * sqrt(5 - 4 * cos(pi / number_of_star_branches)))

% -- To reproduce [BCV2022], Section 6.1.1, geometry \Omega_* --
%  * Case r_* = 1.83e−2, uncomment the line below:
% perimeters = 0.4;
%  * Case r_* = 4.02e−2, uncomment the line below:
perimeters = 0.88;
% --

% -- To study the convergence with respect to the size of the feature --
% perimeters = 1e-2 ./ 2.^(-5:3);
% --

interior_radiuses = perimeters / (2 * number_of_star_branches * sqrt(5 - 4 * cos(pi / number_of_star_branches)));
filename = 'results/test19_neg';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) ones(size(x)); 
problem_data.h = @(x, y, ind) zeros(size(x));

[problem_data, problem_data_0] = set_boundary_conditions(problem_data, number_of_star_branches);
defeaturing_data = set_defeaturing_data(problem_data, number_of_star_branches);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [8 8];
method_data.nquad = [5 5];


%% Main
number_of_interior_radiuses = numel(interior_radiuses);
error_H1s = zeros(1, number_of_interior_radiuses);
error_H1s_boundary_representation = zeros(1, number_of_interior_radiuses);
estimator = zeros(1, number_of_interior_radiuses);
measure_of_gamma = zeros(1, number_of_interior_radiuses);
norm_of_u = zeros(1, number_of_interior_radiuses);
relative_error_H1s = zeros(1, number_of_interior_radiuses);

for iter = 1:number_of_interior_radiuses
    radius_in = interior_radiuses(iter);
    fprintf('----- interior radius = %f -----\n', radius_in);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = build_geometry(radius_in, number_of_star_branches);
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
    save(filename, 'interior_radiuses', 'error_H1s', 'error_H1s_boundary_representation', ...
        'estimator', 'measure_of_gamma', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(interior_radiuses, error_H1s, '+-r', interior_radiuses, estimator, '+-b', ...
           interior_radiuses, interior_radiuses.^2, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^2', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf('For interior_radius = %e, \n', interior_radiuses(end))
fprintf('    * Estimator                                      = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}            = %e \n', error_H1s(end))
fprintf('    * Relative error |u-u_0|_{1,Omega}/|u|_{1,Omega} = %e \n', relative_error_H1s(end))
fprintf('    * Effectivity index                              = %f \n', estimator(end) / error_H1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_F] = build_geometry(radius_in, number_of_star_branches)
    radius_out = 2 * radius_in;
    
    for ii = 0:number_of_star_branches - 1
        theta_in = ii * 2*pi/number_of_star_branches;
        if ii > 0
            line = nrbline([radius_out * cos(theta_out), radius_out * sin(theta_out)], ...
                           [radius_in * cos(theta_in), radius_in * sin(theta_in)]);
            circ = nrbcirc(1, [0,0], theta_out, theta_in);
            srf(2*ii) = nrbruled(line, circ);

            srf_F(ii) = nrbruled(nrbline([radius_in * cos((ii-1) * 2*pi/number_of_star_branches), ...
                                          radius_in * sin((ii-1) * 2*pi/number_of_star_branches)], [0, 0]), line);
        end
        theta_out = theta_in + pi/number_of_star_branches;
        line = nrbline([radius_in * cos(theta_in), radius_in * sin(theta_in)], ...
                       [radius_out * cos(theta_out), radius_out * sin(theta_out)]);
        circ = nrbcirc(1, [0,0], theta_in, theta_out);
        srf(2*ii+1) = nrbruled(line, circ);
    end
    theta_in = 0;
    line = nrbline([radius_out * cos(theta_out), radius_out * sin(theta_out)], ...
                   [radius_in * cos(theta_in), radius_in * sin(theta_in)]);
    circ = nrbcirc(1, [0, 0], theta_out, theta_in);
    srf(2*number_of_star_branches) = nrbruled(line, circ);
    
    srf_F(number_of_star_branches) = nrbruled(nrbline([radius_in * cos((number_of_star_branches-1) * 2*pi/number_of_star_branches), ...
                                                       radius_in * sin((number_of_star_branches-1) * 2*pi/number_of_star_branches)], [0,0]), line);
    for ii = 1:number_of_star_branches
        srf_F(ii) = nrbdegelev(srf_F(ii), [1, 1]);
        srf_F(ii) = nrbkntins(srf_F(ii), {0.5, 0.5});

        srf(2*ii-1) = nrbdegelev(srf(2 * ii - 1), [0, 1]);
        srf(2*ii) = nrbdegelev(srf(2 * ii), [0, 1]);

        extra_knots = linspace(0, 1, 10);
        srf(2*ii-1) = nrbkntins(srf(2 * ii - 1), {0.5, extra_knots});
        srf(2*ii) = nrbkntins(srf(2 * ii), {0.5, extra_knots});
    end

    srf_0 = srf;
    srf_0(2 * number_of_star_branches + (1:number_of_star_branches)) = srf_F;
end

function [problem_data, problem_data_0] = set_boundary_conditions(problem_data, number_of_star_branches)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = 1:2:4 * number_of_star_branches;
    problem_data.drchlt_sides = 2:2:4 * number_of_star_branches;

    % Simplified problem
    problem_data_0.nmnn_sides = []; 
    problem_data_0.drchlt_sides = 1:2 * number_of_star_branches;
end

function defeaturing_data = set_defeaturing_data(problem_data, number_of_star_branches)
    defeaturing_data.omega_star_patches_in_omega0 = 1:2 * number_of_star_branches;

    defeaturing_data.gamma_n_sides_in_omega = 1:2:4 * number_of_star_branches;
    defeaturing_data.gamma_n_sides_in_omega0.patch = 1:2 * number_of_star_branches;
    defeaturing_data.gamma_n_sides_in_omega0.local_side_in_patch = 3 * ones(2 * number_of_star_branches, 1);

    defeaturing_data = decompose_into_single_defeaturing_terms(defeaturing_data, problem_data);
end
