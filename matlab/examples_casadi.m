%%% Examples of Casadi

import casadi.*

select_program = 1;

%%% Wed Nov 22 05:52:55 AM EST 2023
%%% First example from web
if (select_program == 1)

    %% Algorithmic Differentiation (AD)
    % Create scalar/matrix symbols
    x = MX.sym('x',5);

    % Compose into expressions
    y = norm(x,2);

    % Sensitivity of expression -> new expression
    grad_y = gradient(y,x)

    % Create a Function to evaluate expression
    f = Function('f',{x},{grad_y});

    % Evaluate numerically
    grad_y_num = f([1;2;3;4;5]);

    %% Dynamic systems
    x = MX.sym('x',2); % Two states

    % Expression for ODE right-hand side
    z = 1-x(2)^2;
    rhs = [z*x(1)-x(2);x(1)];

    ode = struct;    % ODE declaration
    ode.x   = x;     % states
    ode.ode = rhs;   % right-hand side

    % Construct a Function that integrates over 4s
    F = integrator('F','cvodes',ode,0,4);

    % Start from x=[0;1]
    res = F('x0',[0;1]);

    disp(res.xf)

    % Sensitivity wrt initial state
    res = F('x0',x);
    S = Function('S',{x},{jacobian(res.xf,x)});

    disp(S([0;1]))

    %% Nonlinear and quadratic programming

    % Symbols/expressions
    x = MX.sym('x');
    y = MX.sym('y');
    z = MX.sym('z');
    f = x^2+100*z^2;
    g = z+(1-x)^2-y;

    nlp = struct;            % NLP declaration
    nlp.x = [x;y;z];         % decision vars
    nlp.f = f;               % objective
    nlp.g = g;               % constraints

    % Create solver instance
    F = nlpsol('F','ipopt',nlp);

    % Solve the problem using a guess
    F('x0',[2.5 3.0 0.75],'ubg',0,'lbg',0)

    %% Composition of the above

    x = MX.sym('x',2); % Two states
    p = MX.sym('p');   % Free parameter

    % Expression for ODE right-hand side
    z = 1-x(2)^2;
    rhs = [z*x(1)-x(2)+2*tanh(p);x(1)];

    % ODE declaration with free parameter
    ode = struct('x',x,'p',p,'ode',rhs);

    % Construct a Function that integrates over 1s
    F = integrator('F','cvodes',ode,0,1);

    % Control vector
    u = MX.sym('u',5,1);

    x = [0;1]; % Initial state
    for k=1:5
      % Integrate 1s forward in time:
      % call integrator symbolically
      res = F('x0',x,'p',u(k));
      x = res.xf;
    end

    % NLP declaration
    nlp = struct('x',u,'f',dot(u,u),'g',x);

    % Solve using IPOPT
    solver = nlpsol('solver','ipopt',nlp);
    res = solver('x0',0.2,'lbg',0,'ubg',0);

    plot(full(res.x))
end

%%% Wed Nov 22 05:52:55 AM EST 2023
% 1) simple trajectory optimization
if (select_program == 2)

end