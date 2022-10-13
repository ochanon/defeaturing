clear
close all
clc

eps_values = 1e-2 ./ 2.^(0:6);
filename = 'results/test16_drchlt0';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 10*cos(3*pi*x).*sin(5*pi*y); 
problem_data.h = @(x, y, ind) zeros(size(x));
[problem_data, problem_data_0] = determineBC(problem_data);
ext = 4;

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32];
method_data.nquad = [5 5];

%%
nn = numel(eps_values);
errh1s = zeros(1,nn);
errh1s_interface = zeros(1,nn);
est = zeros(1,nn);
meas_gamma = zeros(1,nn);
normu = zeros(1,nn);
errh1s_rel = zeros(1,nn);

for ii = 1:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    normu(ii) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    errh1s_rel(ii) = errh1s(ii)/normu(ii);
%     [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
%         est_negative(space_0, u_0, problem_data.gamma_sides, problem_data_0.omega_patches, problem_data.g, msh, space, u);

    [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
        est_negative_v2(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end

%%
if saveIt
    save(filename, 'eps_values', 'errh1s', 'errh1s_interface', 'est', 'meas_gamma', 'normu', 'errh1s_rel')
end
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r', eps_values, est, '+-b', eps_values, eps_values, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

%%
function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    s = epsilon/(2*sqrt(2));
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
    srf_0(5) = nrbkntins(srf_0(5), {0.5, (ext-1)*epsilon/(1-epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], 5*pi/4, 7*pi/4);
    crv2 = nrbline([0, 0], [1, 0]);
    srf_0(6) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(6) = nrbkntins(srf_0(6), {0.5, (ext-1)*epsilon/(1-epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], -pi/4, 0);
    crv2 = nrbline([1, 0], [1, 1]);
    srf_0(7) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(7) = nrbkntins(srf_0(7), {0.5, (ext-1)*epsilon/(1-epsilon)});
    
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
%     problem_data.gamma_info = {[2 3 4], {4, 4, 4}, [1 3; 2 3; 3 3]};

    % Simplified problem
    problem_data_0.nmnn_sides = [1 2 3 4 7];
    problem_data_0.drchlt_sides = [5 6 8];
    problem_data_0.omega_patches = 5:7;
    problem_data_0.gamma_sides = cell(7,1); 
    problem_data_0.gamma_sides(5:7) = {3, 3, 3}; % relative to each patch
end