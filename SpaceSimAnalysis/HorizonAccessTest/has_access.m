function tf = has_access(r, t)

% r = system position
% t = target Position
% p = pointing vector

p = t - r;

p_n = norm(p);
r_n = norm(r);
t_n = norm(t);

% cos(alpha) = t*r/t/r
alpha = acos(dot(t*r)/t_n/r_n);
% cos(beta) = p*(-r)/p/r
beta = acos(dot(p, -r)/p_n/r_n);

if (alpha + beta <= pi/2)
    tf = 1;
else
    tf = 0;
end