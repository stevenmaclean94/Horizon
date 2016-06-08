function dy = eoms(t, y)

% The inputs to this funtion are t and y.  t is not used explicitly in
% this function.  The vector r is made up of the position and velocity of
% the spacecraft.  The position, r is made up of members 1, 2 and 3 of the
% vector y.  The velocity, v is made up of members 4, 5, and 6.  The
% velocity is the derivative of the position.

% define mu
mu = 3.98601e5;


% calculate the norm of the position vector, r 
r_norm = norm(y(1:3));

% clear dy
dy = zeros(6,1);    % a column vector

% dy(1:3) are the new velocity terms used in ode45
dy(1) = y(4);
dy(2) = y(5);
dy(3) = y(6);

% dy(4:6) are the new acceleration terms used in ode45
dy(4) = -mu/r_norm^3 * y(1);
dy(5) = -mu/r_norm^3 * y(2);
dy(6) = -mu/r_norm^3 * y(3);