%%%% field_conf.m
clear;
syms t r theta phi g gp rho f A B real

rho = sym('rho(r)')
f = sym('f(r)')
A = sym('A(r)')
B = sym('B(r)')

%% metric
n = [-1,0,0,0;
     0,1,0,0;
     0,0,1,0;
     0,0,0,1]

 %% create array of generators
T = cell(4, 1)
T{1} = [0 1; 1 0]
T{2} = [0 -i; i 0]
T{3} = [1 0; 0 -1]
T{4} = (gp/g)*[1 0; 0 1]

%% unitary transformation
U = i*[cos(theta/2) sin(theta/2)*exp(-i*phi);-sin(theta/2)*exp(i*phi) cos(theta/2)]
Ud = ctranspose(U)

%% higss doublet
xiH = i*[sin(theta/2)*exp(-i*phi);-cos(theta/2)]
xidH = ctranspose(xiH)
xiU = [0;1]
xidU = ctranspose(xiU)
xi = xiU
xid = xidU

%% unit vectors
t_unit = [1;0;0;0]
r_unit = [0;sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)]
theta_unit = [0;cos(theta)*cos(phi);cos(theta)*sin(phi);-sin(theta)]
phi_unit = [0;-sin(phi);cos(phi);0]

%% gauge vector components
A1C = [0; 0; -f*sin(phi)/(g*r); -f*cos(phi)/(g*r)]
A2C = [0; 0; f*cos(phi)/(g*r); -f*sin(phi)/(g*r)]
A3C = [A/g; 0; 0; -(1-cos(theta))/(g*r*sin(theta))]
A4C = [B/gp; 0; 0; -(1-cos(theta))/(gp*r*sin(theta))]

%% gauge vectors
A1 = A1C(1,1)*t_unit + A1C(2,1)*r_unit + A1C(3,1)*theta_unit + A1C(4,1)*phi_unit
A2 = A2C(1,1)*t_unit + A2C(2,1)*r_unit + A2C(3,1)*theta_unit + A2C(4,1)*phi_unit
A3 = A3C(1,1)*t_unit + A3C(2,1)*r_unit + A3C(3,1)*theta_unit + A3C(4,1)*phi_unit
A4 = A4C(1,1)*t_unit + A4C(2,1)*r_unit + A4C(3,1)*theta_unit + A4C(4,1)*phi_unit

%% create array to hold gauge fields
A = cell(4, 1)
A{1} = A1
A{2} = A2
A{3} = A3
A{4} = A4

%% array gauge fields, upper index
Au = cell(4, 1)
Au{1} = n*A1
Au{2} = n*A2
Au{3} = n*A3
Au{4} = n*A4

%% array gauge fields, componets
AC = cell(4, 1)
AC{1} = n*A1C
AC{2} = n*A2C
AC{3} = n*A3C
AC{4} = n*A4C

%% array gauge fields, componets upper index
ACu = cell(4, 1)
ACu{1} = n*A1C
ACu{2} = n*A2C
ACu{3} = n*A3C
ACu{4} = n*A4C

%% create epsilon tensor
epsi = zeros(3, 3, 3)
for j = 1:3
    for k = 1:3
        for l = 1:3
            if j == 1 & k == 2 & l == 3 
                epsi(j, k, l) = 1
            elseif j == 1 & k == 3 & l == 2 
                epsi(j, k, l) = -1
            elseif j == 2 & k == 1 & l == 3 
                epsi(j, k, l) = -1
            elseif j == 2 & k == 3 & l == 1 
                epsi(j, k, l) = 1
            elseif j == 3 & k == 1 & l == 2 
                epsi(j, k, l) = 1
            elseif j == 3 & k == 2 & l == 1 
                epsi(j, k, l) = -1       
            end   
        end
    end
end