% Copyright (C) 2006-2008 Brian C. Fabien. All rights reserved.
%
%  Redistribution and use in source or binary forms, with or without
%  modification, are permitted provided that the following conditions
%  are met:
%  1. Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%  2. Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%  3. All advertising materials mentioning features or use of this software
%     must display the following acknowledgement:
%        This product includes software developed by Brian C. Fabien.
%  4. The name of the author (Brian C. Fabien) may NOT be used to endorse 
%     or promote products derived from this software
%     without specific prior written permission.
% 
%  THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, 
%  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
%  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE 
%  AUTHOR (Brian C. Fabien) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
%  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
%  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
%  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
%  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
%  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% This function solves implicit differential equations (IDE).
% The solution technique is based on the 3-stage Radau IIA implicit Runge-Kutta
% coefficients.
%
% Usage:
%
%    [TOUT,YOUT,INFO] = ride('IDEFUN', 'IDEJAC', TSPAN, Y0, YP0)  or
%    [TOUT,YOUT,INFO] = ride('IDEFUN', 'IDEJAC', TSPAN, Y0, YP0, OPTIONS) or
%    [TOUT,YOUT,INFO] = ride('IDEFUN', '', TSPAN, Y0, YP0)  or
%    [TOUT,YOUT,INFO] = ride('IDEFUN', '', TSPAN, Y0, YP0, OPTIONS)
%
% Here,
%
% TOUT     is a (Ntspan x 1) vector of output times at which the solution is 
%          computed. (Ntspan is defined below.)
%
% YOUT     is a (Ntspan x n) matrix of solutions to the IDE, where n is the
%          dimension of the IDE.  (Ntspan is defined below.)
%
% INFO     is a structure that provides some statistics related to the solution
%          algorithm.  In particular,
%
%               INFO.nfun    - is the number of function evaluations.
%               INFO.njac    - is the number of Jacobian evaluations.
%               INFO.naccept - is the number of successful steps.
%               INFO.nreject - is the number of failed steps.
%
% IDEFUN   is a function that returns the (n x 1) vector res = phi(y, yp, t), 
%          where y is the (n x 1) state vector of the IDE, yp = dy/dt is the 
%          (n x 1) state derivative, t is the time, and phi is the (n x 1)
%          system of IDE. This function has the form
%          
%               function res = IDEFUN(y, yp, t)
%
% IDEJAC   is a function that computes the (n x n) Jacobian J = phi/dy, and the
%          (nxn) matrix M = dphi/dyp.  This function has the form
%                  
%               function [J, M] = IDEJAC(y, yp, t)
%
%          where
%
%                |dphi_1/dy_1 dphi_1/dy_2 ... dphi_1/dy_n|
%                |dphi_2/dy_1 dphi_2/dy_2 ... dphi_2/dy_n|
%           J =  |   .                                   | and
%                |   .                                   |
%                |dphi_n/dy_1 dphi_n/dy_2 ... dphi_n/dy_n|
%
%                |dphi_1/dyp_1 dphi_1/dyp_2 ... dphi_1/dyp_n|
%                |dphi_2/dyp_1 dphi_2/dyp_2 ... dphi_2/dyp_n|
%           M =  |   .                                      |
%                |   .                                      |
%                |dphi_n/dyp_1 dphi_n/dyp_2 ... dphi_n/dyp_n|
%
%          If IDEJAC is not included in the argument list for ride then J and M
%          are computed via a finite difference approximation.
%
% TSPAN    This vector defines the limits of integration, as well as 
%          intermidiate time points where the solution to the IDE is computed.
%          Specifically,
%
%          TSPAN is a (Ntspan x 1) vector
%
%               TSPAN = [t_1; t_2; t_3; ...; t_Ntspan]
%
%          where the times t_1, t_2, ..., t_Ntspan is a monotone increasing 
%          (or decreasing) sequence.  The limits of integration are t_1 (the
%          initial time, and T_Ntspan (the final time).  Solutions to the IDE
%          are computed a the points t_2, t_3, ..., t_Ntspan, and stored in the
%          matrix YOUT.  On successful termination of the function TOUT = TSPAN.
%
% Y0       This (n x 1) vector gives the initial state for the IDE.
%
% YP0      This (n x 1) vector gives the initial state derivative for the IDE.
%          Moreover, it is assumed that these initial conditions are consistent,
%          i.e., res = phi(Y0, YP0, TSPAN(1)) = 0.  Note that the function 
%          terminates if ||phi(Y0, YP0, TSPAN(1))|| > sqrt(eps), where eps is 
%          the machine precision.
%
% OPTIONS  This is a structure that provides various options for the numerical
%          solution algorithm.  The default values for the members of this 
%          structure can be ascertained using the function 
%
%                OPTIONS = ide_options()
%
%          The function options can be assigned using the function 
%
%                OPTIONS = ide_set_option(OPTIONS, 'OPTION', value)
%
%          The options for the function ride and the default values are as 
%          follows;
%
%          OPTIONS.ATOL 
%             is the absolute error tolerance. (Default 1.0e-6).
%
%          OPTIONS.RTOL 
%             is the relative error tolerance. (Default 1.0e-6).
%
%          Note that ATOL and RTOL can also be (n x 1) vectors.
%
%          OPTIONS.INITIAL_STEP_SIZE 
%             is the initial step size. (Default 1.0e-8).
%
%          OPTIONS.MAX_STEPS 
%             is maximum number of steps that is performed by the function. 
%             (Default 1000).
%
%          OPTIONS.DIFF_INDEX 
%             is a (n x 1) vector that indicates the differentiation index of
%             the i-th state variable, y(i). (Default [], i.e., all state 
%             variables are assumed to be index 0).
%
%          OPTIONS.T_EVENT 
%             is a (P x 1) vector of time points of the form
%
%               OPTIONS.T_EVENT = [te_1;te_1;te_2;...;te_P]
%
%          where the times te_1, te_2, ..., te_P is a monotone increasing 
%          (or decreasing) sequence.  In the algorithm the step sizes are 
%          selected such that the event times, T_EVENT(k), are at the end 
%          (or start) of the integration interval.  This is useful for modeling
%          discontinuous time varying inputs or, discontinuous implicit 
%          differential equations.
%
% Example:
%
% (1) The file ide.m contains
%
% function res = ide(y, yp, t)
%
% res = zeros(2, 1);
%
% if ((t >= 0) && (t < 1))
%    u = 1;
% elseif ((t >= 1) && (t < 2))
%    u = -1;
% elseif((t >= 2) && (t < 3))
%    u = 1;
% else
%    u = 0;
% end;
% res(1) = yp(1) - y(2);
% res(2) = yp(2) + 0.1 * y(2) + 3.0 * y(1) - u;
% return;
%
% (2) The file idejac.m contains
%
% function [J, M] = idejac(y, yp, t)
%
% J = [0, -1; 0.1, 3.0];
% M = eye(2);
% return;
%
% (3) The file runide.m contains
%
% clear all;
% tspan = linspace(0, 60.0, 600)';
% ride_opt = ide_options();
% ride_opt = ide_set_option(ride_opt, 'T_EVENT', [1; 2; 3]);
% y0 = [0; 0];
% yp0 = [0; 0];
% yp0 = -ide(y0, yp0, 0);
%
% [T, Y, info] = ride('ide', 'idejac', tspan, y0, yp0, ride_opt);
%
%
% See also, ide_options, ide_set_option
% 

function [Tout, Yout, info] = ride(IDEFUN, IDEJAC, tspan, y0, yp0, options)

global ride_data;

if (nargin < 6)
    RTOL = 1.0e-6;
    ATOL = 1.0e-6;
    MAX_STEPS = 1000;
    DIFF_INDEX = [];
    INITIAL_STEP_SIZE = 1.0e-8;
    T_EVENT = [];
else
    RTOL = options.RTOL;
    ATOL = options.ATOL;
    MAX_STEPS = options.MAX_STEPS;
    DIFF_INDEX = options.DIFF_INDEX;
    INITIAL_STEP_SIZE = options.INITIAL_STEP_SIZE;
    T_EVENT = options.T_EVENT;
end;

% The Runge-Kutta coefficients: the 3-stage Radau IIa method

s = 3; % The number of stages

A = [0.1968154772236604, -0.0655354258501984,  0.0237709743482202; 
     0.3944243147390873,  0.2920734116652284, -0.0415487521259979; 
     0.3764030627004673,  0.5124858261884216,  0.1111111111111111];
          
%b = A(3,:)';

c = [0.155051025721682;
     0.644948974278318;
     1.0];

% Transformation of inv(A) to Jordan canonical form

% inv(T) * inv(A) * T = [mu1  0    0;
%                        0   mu2  mu3;
%                        0  -mu3  mu2];
                  
% Ai = inv(A);

Ai = [3.224744871391587,  1.167840084690408, -0.253197264742182; 
     -3.567840084690405,  0.775255128608409,  1.053197264742182; 
      5.531972647421809, -7.531972647421808,  4.999999999999999];

T = [0.091232394870893, -0.128458062178301,  0.027308654751321; 
     0.241717932707107,  0.185635951030957, -0.348248904396575; 
     0.966048182615093,  0.909403517646862,  0.0];

% Ti = inv(T);    
Ti = [4.325579890063155,  0.339199251815810, 0.541770539935875; 
     -4.595010367196070, -0.360327197335856, 0.524105686036759; 
      0.552969749058155, -2.828147131551269, 0.655417747196001];
      
mu = [3.63783425274449; 
      2.68108287362775;
      3.05043019924741];

mu23 = mu(2) - sqrt(-1) * mu(3);
  
c2mc1 = c(2) - c(1);
c3mc1 = c(3) - c(1);
c3mc2 = c(3) - c(2);

% Coefficients used to estimate the local truncation error
% J. de Swart and G. Soderlind, On the construction of error estimators for 
% implicit Runge-Kutta methods, Journal of Computational and Applied 
% Mathematics, 87, 347--358, 1997.

b0 = 0.02;
v = [0.0311615640944986; -0.0178282307611655; 0.2815554962623449];

MAX_NEWT_ITER = 10; % The maximum number of Newton iterations per step.

% This data structure is used to store the problem information

ride_data = struct('A', A, 'Ai', Ai, 'c', c, 'mu23', mu23, ...
    'c2mc1', c2mc1, 'c3mc1', c3mc1, 'c3mc2', c3mc2, 'D', [], ...
    'T', T, 'Ti', Ti, 'v', v, 'b0', b0, 'n', 0, 'sn', 1, 's3n', 1, ...
    'rtol', RTOL, 'atol', ATOL, 'NTOL', 0.01, ...
    'MAX_NEWT_ITER', MAX_NEWT_ITER, 'M', [], 'J', [], ...
    'L', [], 'U', [], 'P', [], 'nnewton', 0, ...
    'nfun', 0, 'njac', 0, 'n_newt', 0, 'nr_last', 1, 'h_last', 0, ...
    'DIFF_INDEX', DIFF_INDEX, 'error_scale_index', [], 'FD_JACOBIAN', 0, ...
    'Z',[], 'y0', [], 'ym1', [], 'mu', mu, ...
    'ifail', 0, 'ifail_count', 0, 'no_update', 0, 'first_step', 1);

% Some statistical information is stored in the info data structure
    
info = struct('nfun', 0, 'njac', 0, 'naccept', 0, 'nreject', 0, 'nnewton', 0);

Tout = []; % The output time points
Yout = []; % The output states

% Check the input data

[n, m] = size(y0);
if (m > 1)
    printf('Error: expecting an n by 1 matrix for y0\n');
    return;
end;

[np, mp] = size(yp0);
if ((np ~= n) || (mp > 1))
    printf('Error: expecting an n by 1 matrix for yp0\n');
    return;
end;

ride_data.n = n;

if (isempty(ride_data.DIFF_INDEX))
    ride_data.DIFF_INDEX = zeros(n, 1);
else
    [nd, md] = size(ride_data.DIFF_INDEX);
    if ((nd ~= n) || (md > 1))
        printf('Error: expecting an n by 1 matrix for DIFF_INDEX\n');
        return;
    end;
    if (min(ride_data.DIFF_INDEX) < 0)
        printf('Error: invalid data in the matrix DIFF_INDEX\n')
        return;
    end;
end;

index = max(ride_data.DIFF_INDEX - 1, 0);

if ((length(ride_data.rtol) ~= 1) ...
    && ((length(ride_data.rtol) ~= length(ride_data.atol)) ...
        || (length(ride_data.rtol) ~= ride_data.n)))
    printf('Error: invalid dimensions for RTOL and ATOL\n')
    return;
end;

if (min(ride_data.rtol) <= 0)
    printf('Error: invalid input for RTOL\n')
    return;
end;

if (min(ride_data.atol) <= 0)
    printf('Error: invalid input for ATOL\n')
    return;
end;

if (strcmp(IDEJAC,'') == 1)
    ride_data.FD_JACOBIAN = 1;
end;

% Is tspan in the correct order?

ntspan = length(tspan);

if (ntspan < 2)
    printf('Error: TSPAN should be at least a 2 x 1 vector\n')
    return;
end;

t0 = tspan(1);
tend = tspan(ntspan);
dt = tend - t0;
if (abs(dt) <= (10.0 * eps))
    printf('Error: |t_final - t_initial| too small\n');
    return;
end;

if (min(dt * (tspan(2 : ntspan) - tspan(1 : (ntspan - 1)))) <= 0)
    printf('Error: TSPAN not monotone increasing (or monotone decreasing)\n');
    return;
end;

h = min(abs(dt), abs(INITIAL_STEP_SIZE));
inext = 2;
tnext = tspan(inext);

% Are the time events in the correct order?

ntevent = length(T_EVENT);
next_event = 1;

if (ntevent > 1)
    if (min(dt * (T_EVENT(2 : ntevent) - T_EVENT(1 : (ntevent - 1)))) <= 0)
        printf('Error: T_EVENT and TSPAN must be both monotone increasing (or both monotone decreasing)\n');
        return;
    end;    
end;

if (ntevent > 0)
    dt = T_EVENT(next_event) - t0;
    if (abs(dt) <= 10.0 * eps)
        printf('Error: T_EVENT(1) is too close to t0');
        return;
    end;
    h = min(abs(dt), h);
end;

% Check for backwards integration

if ((h * (tend - t0)) < 0) 
    h = -h;
end;

% Set the Newton iteration convergence tolerance

ctol = 0.01;
if (length(ride_data.rtol) == 1)
    ride_data.NTOL = max(10.0 * eps / ride_data.rtol, ...
                         min(ctol, sqrt(ride_data.rtol)));
else
    rtol = min(ride_data.rtol);
    ride_data.NTOL = max(10.0 * eps / rtol, min(ctol, sqrt(rtol)));
end;
% For very strict convergence of the Newton iteration set
% ride_data.NTOL very small, say <= 0.001

% Set the error norm scaling factor

error_scale =  abs(y0) .* ride_data.rtol + ride_data.atol;
ride_data.error_scale_index = error_scale ./ (h .^ index);

% Are the initial conditions consistent?

phi0 = feval(IDEFUN, y0, yp0, t0);
ride_data.nfun = ride_data.nfun + 1;
norm_phi0 = norm(phi0, 'fro') / ride_data.sn;
if (norm_phi0 > sqrt(eps))
    %printf('Error: the initial conditions are not consistent: |phi0| = %e\n', norm_phi0);
    %return;
end;

% Compute the Jacobian at t0

if (ride_data.FD_JACOBIAN == 0)
    [ride_data.J, ride_data.M] = feval(IDEJAC, y0, yp0, t0);
else
    [ride_data.J, ride_data.M] = ride_fd_jacobian(IDEFUN, y0, yp0, t0);
    ride_data.nfun = ride_data.nfun + 2 * n + 1;
end;
ride_data.njac = ride_data.njac + 1;

% Initialize some parameters in the data structure

t = t0;
y = y0;
yp = yp0;

ride_data.sn = sqrt(n);
ride_data.s3n = sqrt(3.0 * n);
ride_data.h_last = h;
ride_data.y0 = y0;
ride_data.ym1 = y0;
ride_data.Z = zeros(n, s);

Tout = zeros(ntspan, 1);
Yout = zeros(ntspan, n);
Tout(1) = t0;
Yout(1, :) = y';
Yp_last = [];

for step = 1 : MAX_STEPS

    % Newton iteration
    if (isempty(Yp_last))
        [Yp0, Z0] = formZ(h, h, [yp0, yp0, yp0]);
    else
        [Yp0, Z0] = formZ(h, ride_data.h_last, Yp_last);
    end;

    [Y, Yp, newton_iteration_error] = ride_newton(IDEFUN, t, h, y, Yp0, Z0);
    
    if (newton_iteration_error == 0) % The Newton iteration converged
        
        % Compute the local truncation error
        [hnew, lte_error] = ride_lte(IDEFUN, Y, Yp, yp, t, h);
        
        if (lte_error == 0) % The local error is acceptable
        
            info.naccept = info.naccept + 1;
            y1  = Y(:, 3);
            yp1 = Yp(:, 3);
            Yp_last = Yp;
            t1  = t + h;
            
            % compute the coefficients for the interpolation formula
            ride_data.D = ride_newton_coefficients(ride_data.Z);
            
            % Perform dense output via a Newton interpolation formula 
            % if tnext is in the interval t < tnext <= t1, 
            % (OR t > tnext >= t1 if we are integrating backwards)
            tau = (tnext - t0) / h;
            while (tau <= 1.0)
                ynext = ride_newton_interp(tau, y);
                Tout(inext) = tnext;
                Yout(inext, :) = ynext';
                inext = inext + 1;
                if (inext <= ntspan) 
                    tnext = tspan(inext); 
                else
                    tnext = tend + 10.0 * h;  % put tnext past the final time
                end;
                tau = (tnext - t0) / h;
            end;
            
            if (abs(tend - t1) < 10.0 * eps) % t1 is at tend
                info.nfun = ride_data.nfun;
                info.njac = ride_data.njac;
                info.nnewton = ride_data.nnewton;
                return;
            end;
            
            if (ride_data.first_step == 1)
                ride_data.first_step = 0;
            end;
            ride_data.y0 = y1;
            ride_data.ym1 = y;
            y = y1;
            yp = yp1;
            t0 = t1;
            t = t0;
            h = hnew;
            
            if ((tend - (t + h)) * h < 0)
                h = tend - t;
            end;

            if (next_event <= ntevent)
                dt = T_EVENT(next_event) - t;
                if (abs(dt) < 10.0 * eps)
                    next_event = next_event + 1;
                    if (next_event > ntevent)
                        t_event_n = tend + 10.0 * h;
                    else
                        t_event_n = T_EVENT(next_event);
                    end;
                else
                    t_event_n = T_EVENT(next_event);                    
                end;
                
                % Do not integrate past the next event
                if ((t_event_n - (t + h)) * h < 10.0 * eps)
                    h = t_event_n - t;
                end;
            end;
            if (ride_data.no_update == 0)
                if (ride_data.FD_JACOBIAN == 0)
                    [ride_data.J, ride_data.M] = feval(IDEJAC, y, yp, t);
                else
                    [ride_data.J, ride_data.M] = ride_fd_jacobian(IDEFUN, y, yp, t);
                    ride_data.nfun = ride_data.nfun + 2 * n + 1;
                end;
                ride_data.njac = ride_data.njac + 1;
            end;
            
            % Keep the step size the same for at least three steps
            % after a failure of the Newton iteration, or an unsatisfactory LTE
            if (ride_data.ifail == 1) 
                if (ride_data.ifail_count == 3)
                    ride_data.ifail = 0;
                else
                    ride_data.ifail_count = ride_data.ifail_count + 1;
                end;
            end;
            
        else
            info.nreject = info.nreject + 1;
            h = hnew;
            if (abs(h) <= 10.0 * eps)
                printf('ride: step size too small: 0\n');
                info.nfun = ride_data.nfun;
                info.njac = ride_data.njac;
                info.nnewton = ride_data.nnewton;
                return;
            end;
            ride_data.ifail = 1;
            ride_data.ifail_count = 1;
        end;
    else
        info.nreject = info.nreject + 1;
        h = h * 0.5;
        if (abs(h) <= 10.0 * eps)
            printf('ride: step size too small: 1\n');
            info.nfun = ride_data.nfun;
            info.njac = ride_data.njac;
            info.nnewton = ride_data.nnewton;
            return;
        end;
        ride_data.ifail = 1;
        ride_data.ifail_count = 1;
    end;

    % set the error norm scaling factor
    error_scale =  abs(y) .* ride_data.rtol + ride_data.atol;
    ride_data.error_scale_index = error_scale ./ (h .^ index);
    
end;
printf('ride: too many steps\n');
info.nfun = ride_data.nfun;
info.njac = ride_data.njac;
info.nnewton = ride_data.nnewton;
return;

function [hnew, err] = ride_lte(IDEFUN, Y, Yp, yp0, t0, h)
% This function computes the local truncation error and a new step size

global ride_data;

beta = 0.9;
fac0 = 0.2;
fac1 = 5;
err = 0;

y1 = Y(:, 3);
t1 = t0 + h;

Yp1 = Yp * ride_data.v;

yp1 = ride_data.mu(1) * (Yp1 - ride_data.b0 * yp0);

g = feval(IDEFUN, y1, yp1, t1);
ride_data.nfun = ride_data.nfun + 1;

r = -(ride_data.U \ (ride_data.L \ (ride_data.P * g)));

nr = norm(r ./ ride_data.error_scale_index, 'fro') / ride_data.sn;

beta1 = beta * (2.0 * ride_data.MAX_NEWT_ITER + 1.0) ...
        / (2.0 * ride_data.MAX_NEWT_ITER + ride_data.n_newt);
fac2 = min(fac1, max(fac0, beta1 * nr ^ (-0.25)));
if (nr <= eps)
    fac3 = fac1;
else
    nr1 = ((1.0 / nr) ^ 0.25) * (h / ride_data.h_last) ...
     * ((ride_data.nr_last / nr) ^ 0.25);
    fac3 = min(fac1, max(fac0, beta1 * nr1));
end;

fac = min(fac2, fac3);

if ((ride_data.n_newt == 1) && (fac >= 0.9) && (fac <= 1.1))
    fac = 1.0;
    ride_data.no_update = 1;
else
    ride.no_update = 0;
end;

if ((ride_data.ifail == 1) && (fac > 1.0))
    fac = 1;
end;

if (nr >= 1)
    err = -1;
else
    ride_data.h_last = h;
    ride_data.nr_last = nr;
end;

hnew = h * fac;
return;

function [J, M] = ride_fd_jacobian(IDEFUN, y0, yp0, t0)
% This function computes the Jacobian of the IDE using a forward difference approximation

global ride_data;

DBL_EPSILON = eps;
SQRT_DBL_EPSILON = sqrt(eps);

n = ride_data.n;
y1 = y0;
yp1 = yp0;
phi0 = feval(IDEFUN, y0, yp0, t0);
M = zeros(n, n);
J = zeros(n, n);
for i = 1 : n
    h = abs(y0(i)) * SQRT_DBL_EPSILON; 
    if (h <= DBL_EPSILON)
        h = SQRT_DBL_EPSILON;
    end;
    y1(i) = y1(i) + h;
    h = y1(i) - y0(i);
    phi1 = feval(IDEFUN, y1, yp0, t0);
    J(1 : n, i) = ((phi1 - phi0) ./ h);
    y1(i) = y0(i);

    h = abs(yp0(i)) * SQRT_DBL_EPSILON; 
    if (h <= DBL_EPSILON)
        h = SQRT_DBL_EPSILON;
    end;
    yp1(i) = yp1(i) + h;
    h = yp1(i) - yp0(i);
    phi1 = feval(IDEFUN, y0, yp1, t0);
    M(1 : n, i) = ((phi1 - phi0) ./ h);
    yp1(i) = yp0(i);
end;

return;

function D = ride_newton_coefficients(Z)
% Compute the coefficients used in the interpolation formula

global ride_data;

D = zeros(ride_data.n, 3);

D(:, 1) = (1.0 / ride_data.c(1)) * Z(:, 1);
D(:, 2) = (1.0 / ride_data.c2mc1) * ((1.0 / ride_data.c(2)) * Z(:,2) - D(:, 1));
D(:, 3) = (1.0 / (ride_data.c3mc1 * ride_data.c3mc2)) * (Z(:, 3) - D(:, 1) - D(:, 2) * (1.0 - ride_data.c(1)));

return;

function ynext = ride_newton_interp(tau, y)
% This function computes y at tau using a Newton interpolation formula

global ride_data;

ynext = tau * (ride_data.D(:, 1) + ...
            (tau - ride_data.c(1)) * (ride_data.D(:, 2) + ...
            (tau - ride_data.c(2)) * ride_data.D(:, 3)));
ynext = ynext + y;

return;

function [Yp, Z] = formZ(h, h0, Yp0)

global ride_data;
    
t = h / h0;

E = zeros(3, 3);
E(1, 1) = (1/36)*t*(2/5-(1/10)*sqrt(6))*(-6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(2/5-(1/10)*sqrt(6))^2*(-6+sqrt(6))*sqrt(6);
E(1, 2) = -(1/36)*t*(2/5-(1/10)*sqrt(6))*(6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(2/5-(1/10)*sqrt(6))^2*(6+sqrt(6))*sqrt(6);
E(1, 3) = 1+4*t*(2/5-(1/10)*sqrt(6))+(10/3)*t^2*(2/5-(1/10)*sqrt(6))^2;
E(2, 1) = (1/36)*t*(2/5+(1/10)*sqrt(6))*(-6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(2/5+(1/10)*sqrt(6))^2*(-6+sqrt(6))*sqrt(6);
E(2, 2) = -(1/36)*t*(2/5+(1/10)*sqrt(6))*(6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(2/5+(1/10)*sqrt(6))^2*(6+sqrt(6))*sqrt(6);
E(2, 3) = 1+4*t*(2/5+(1/10)*sqrt(6))+(10/3)*t^2*(2/5+(1/10)*sqrt(6))^2;
E(3, 1) = (1/36)*t*(-6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(-6+sqrt(6))*sqrt(6);
E(3, 2) = -(1/36)*t*(6+sqrt(6))^2*sqrt(6)-(5/18)*t^2*(6+sqrt(6))*sqrt(6);
E(3, 3) = 1+4*t+(10/3)*t^2;

Yp = (E * Yp0')';
Z = h * (ride_data.A * Yp')';

return;


function [Y, Yp, err] = ...
    ride_newton(IDEFUN, t, h, y, Ype, Z1e)
% This function performs the Newtons iterations
% Input
%    IDEFUN: The function that evaluates implict differential equations
%    t:     The current time point
%    h:     The increment in t
%    y:     The state at time t
%    Ype:   The state derivatives at time t
%    Z1e:   An estimate of the increments in the stage values
%
% Output
%    Y:     The stage values
%    Yp:    The stage derivatives
%    Z:     The increment in the stage values
%    iter:  The number of iterations performed
%    err:   An error flag: 
%             0 = normal completion
%            -1 = divergence
%            -2 = convergence too slow
%            -3 = too many iterations
%            -4 = singular coefficient matrix

global ride_data;

norm0 = 1;
n = ride_data.n;

H = (ride_data.mu(1) / h) * ride_data.M + ride_data.J;
[L, U, P] = lu(H);
if ((min(abs(diag(U))) <= 10.0 * eps) || (min(abs(diag(L))) <= 10.0 * eps))
    printf('singular H: t = %e\n', t);
    Y = zeros(n, 3);
    Yp = zeros(n, 3);
    ride_data.Z = zeros(n, 3);
    ride_data.n_newt = 0;
    err = -4;
    return;
end;

ride_data.L = L;
ride_data.U = U;
ride_data.P = P;

h1 = (1.0 / h);
H1 = h1 * ride_data.mu23 * ride_data.M + ride_data.J;
[L1, U1, P1] = lu(H1);
if (min(abs(diag(U1))) < 20.0 * eps)
    printf('singular H1: t = %e\n',t);
    Y = zeros(n, 3);
    Yp = zeros(n, 3);
    ride_data.Z = zeros(n, 3);
    ride_data.n_newt = 0;
    err = -4;
    return;
end;

R = zeros(n, 3);
v = zeros(n, 3);

tc1 = t + ride_data.c(1) * h;
tc2 = t + ride_data.c(2) * h;
tc3 = t + ride_data.c(3) * h;

Yp = Ype;
Z  = Z1e;
Y0 = [y, y, y];
Y  = Y0 + Z;

for iter = 1 : ride_data.MAX_NEWT_ITER
    ride_data.nnewton = ride_data.nnewton + 1;
    R(:, 1) = feval(IDEFUN, Y(:, 1), Yp(:, 1), tc1);
    R(:, 2) = feval(IDEFUN, Y(:, 2), Yp(:, 2), tc2);
    R(:, 3) = feval(IDEFUN, Y(:, 3), Yp(:, 3), tc3);
    ride_data.nfun = ride_data.nfun + 3;
    
    r = (ride_data.Ti * R')';
    
    v(:, 1) = -(U \ (L \ (P * r(:, 1))));
    w23 = -(U1 \ (L1 \ (P1 * (r(:, 2) + sqrt(-1) * r(:, 3)))));
    v(:, 2) = real(w23);
    v(:, 3) = imag(w23);
    
    dZ = (ride_data.T * v')';   
    Z = Z + dZ;    
    Y = Y0 + Z; 
    Yp = h1 * (ride_data.Ai * Z')';
    
    e1 = norm([dZ(:, 1); dZ(:, 2); dZ(:, 3);] ./ ...
            [ride_data.error_scale_index; ...
             ride_data.error_scale_index; ...
             ride_data.error_scale_index], 'fro') / ride_data.s3n;

    if (iter == 1)
        if (e1 <= ride_data.NTOL)
            err = 0;
            ride_data.Z = Z;
            ride_data.n_newt = iter;
            return;
        end;
        norm0 = e1;
    else
        theta = e1 / norm0;        
        if (theta >= 1.0) % Newton iteration is diverging
            err = -1;
            ride_data.Z = Z;
            ride_data.n_newt = iter;
            return;
        end;
        
        eta = theta / (1.0 - theta);
        rs = eta * e1;
        if ((rs <= ride_data.NTOL) || (e1 <= 100.0 * eps))
            err = 0;
            ride_data.Z = Z;
            ride_data.n_newt = iter;
            return;
        end;
        
        rs = theta ^ (ride_data.MAX_NEWT_ITER - iter);
        rs = rs / (1.0 - theta);
        rs = rs * e1;
        if (rs > ride_data.NTOL) % Convergence rate is too slow
            ride_data.Z = Z;
            ride_data.n_newt = iter;
            err = -2;
            return;
        end;
        norm0 = e1;        
    end;
end;

err = -3;

return;
     