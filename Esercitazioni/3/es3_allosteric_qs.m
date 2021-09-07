clear all
close all

syms x
syms y
syms z
syms e0
syms k1
syms km1
syms k3
syms km3
syms k2
syms s
syms i

eqn = [k1*e0*s - (k1*s+km1+k2+k3*i)*x - k1*s*y - (k1*s-km3)*z == 0;
    k3*e0*i-k3*i*x - (k3*i+km3+k1*s)*y - (k3*i-km1)*z == 0;
    k3*i*x + k1*s*y - (km1+km3)*z == 0];

sol = solve(eqn,[x y z]);
x_qs = simplify(sol.x)
y_qs = simplify(sol.y)
z_qs = simplify(sol.z)
Vel_qs = k2*x_qs