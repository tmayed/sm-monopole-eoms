%%%% stress_energy.m
field_conf;

%% Calculate abelian field tensor
F_lower_abel = sym('F_lower_abel',[4 4])
for i = 1:4
    for j = 1:4
        F_lower_abel(i, j) = simplify(vect4("grad_i", A{4}(j), i) - vect4("grad_i", A{4}(i), j))
    end
end

F_upper_abel = sym('F_upper_abel',[4 4])
F_half_abel = sym('F_upper_abel',[4 4])
for i = 1:4
    for j = 1:4
        F_upper_abel(i, j) = n(i,i)*n(j,j)*F_lower_abel(i,j)
        F_half_abel(i, j) = n(i,i)*F_lower_abel(i,j)
    end
end

%% Calculate non abelian field tensor
F_lower_non_abel = sym('F_lower_non_abel',[4 4 3])
for i = 1:4
    for j = 1:4
        for k = 1:3
            hold = 0
            for l = 1:3
                for m = 1:3
                    hold = simplify(hold + g*epsi(k, l, m)*(A{l}(i))*(A{m}(j)))
                end
            end
            F_lower_non_abel(i, j, k) = simplify(vect4("grad_i", A{k}(j), i) - vect4("grad_i", A{k}(i), j) + hold)
        end
    end
end

F_upper_non_abel = sym('F_upper',[4 4 3])
F_half_non_abel = sym('F_half',[4 4 3])
for i = 1:4
    for j = 1:4
        for k = 1:3
            F_upper_non_abel(i,j,k) = n(i,i)*n(j,j)*F_lower_non_abel(i,j,k)
            F_half_non_abel(i,j,k) = n(i,i)*F_lower_non_abel(i,j,k)
        end
    end
end

%% contraction of langranian rho term
hold = 0
for i = 1:4
    hold = hold + simplify(vect4("grad_i", rho, i)*vect4("grad_i", rho, i, 1))
end
rho_cont = simplify(hold)

%% contraction of langranian xi term
hold = 0
for i = 1:4
    for j = 1:4
        for k = 1:4
            hold = hold + simplify(A{j}(i)*Au{k}(i)*xid*T{j}*T{k}*xi)
        end
    end
end
xi_cont = simplify(rho^2*hold)

%% contraction of langranian abelian field tensor
F_abel_cont = 0
for i = 1:4
    for j = 1:4
        F_abel_cont =  F_abel_cont + simplify(F_lower_abel(i, j)*F_upper_abel(i, j))
    end
end
F_abel_cont = simplify(-(1/4)*F_abel_cont)

%% contraction of langranian non-abelian field tensor
F_non_abel_cont = 0
for i = 1:4
    for j = 1:4
        for k = 1:3
            F_non_abel_cont =  F_non_abel_cont + simplify(F_lower_non_abel(i, j, k)*F_upper_non_abel(i, j, k))
        end
    end
end
F_non_abel_cont = simplify(-(1/4)*F_non_abel_cont)

%% put all together
lang = simplify(rho_cont + xi_cont + F_abel_cont + F_non_abel_cont)

%% contract rho expression
T_rho = sym('T_rho',[4 4])
for i = 1:4
    for j = 1:4
        T_rho(i, j) = simplify(vect4("grad_i", rho, i, 1)*vect4("grad_i", rho, j, 1))
    end
end
T_rho_ene = T_rho(1, 1)

%% contract xi expression
T_xi= sym('T_xi',[4 4])
for i = 1:4
    for j = 1:4
        for k = 1:4
            for l = 1:4
                T_xi(i, j) = rho^2*simplify(Au{k}(i)*Au{l}(j)*xid*T{k}*T{l}*xi)
            end
        end
    end
end
T_xi_ene = T_xi(1, 1)

%% contract abel expression
T_abel = sym('T_abel',[4 4])
for i = 1:4
    for j = 1:4
        hold = 0
        for k = 1:4
            hold = hold + simplify(F_upper_abel(i,k)*F_half_abel(j,k))
        end
        T_abel(i, j) = simplify(hold)
    end
end
T_abel_ene = T_abel(1, 1)

%% contract non-abel expression
T_non_abel = sym('T_non_abel',[4 4])
for i = 1:4
    for j = 1:4
        hold = 0
        for k = 1:3
            for l = 1:4
                hold = hold + simplify(F_upper_non_abel(i,l,k)*F_half_non_abel(j,l,k))
            end
        end
        T_non_abel(i, j) = simplify(hold)
    end
end
T_non_abel_ene = T_non_abel(1, 1)

%% Whole Stress Energy Tensor
T_StressEnergy = sym('T_StressEnergy',[4 4])
for i = 1:4
    for j = 1:4
        T_StressEnergy(i, j) = simplify(T_rho(i, j) + T_xi(i, j) + T_abel(i, j) + T_non_abel(i, j) + n(i, j)*lang)
    end
end
T_Energy = T_StressEnergy(1, 1)