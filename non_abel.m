%%%% non_abel.m
field_conf;

%% Calculate the 4-derivative of non abelian field tensor
D_F_upper_non_abel = sym('D_F_upper_non_abel',[4 3])
for a = 1:3
    lap = vect4("lap4", Au{a})
    div_grad = vect4("grad", vect4("div", ACu{a}), 1) 
    for v = 1:4
        hold = 0
        for j = 1:4
            for b = 1:3
                for c = 1:3
                    if epsi(a, b, c) ~= 0
                        hold = hold + g*epsi(a, b, c)*(vect4("grad_i", (Au{b}(j))*(Au{c}(v)), j))
                    end
                end 
            end
        end
        D_F_upper_non_abel(v, a) = simplify(lap(v) - div_grad(v) + hold)
    end    
 end 

%% Calculate non abelian field tensor
F_upper_non_abel = sym('F_upper_non_abel',[4 4 3])
for v = 1:4
    for j = 1:4
        for a = 1:3
            hold = 0
            for b = 1:3
                for c = 1:3
                    if epsi(a, b, c) ~= 0 
                        hold = simplify(hold + g*epsi(a, b, c)*(Au{b}(v))*(Au{c}(j)))
                    end
                end
            end
            F_upper_non_abel(v, j, a) = simplify(vect4("grad_i", Au{a}(j), v, 1) - vect4("grad_i", Au{a}(v), j, 1) + hold)
        end
    end
end

%% Calculate the epsilon contractions with the non abelian field tensor
epsi_F_upper_non_abel = sym('epsi_F_upper_non_abel',[4 3])
for a = 1:3
    for v = 1:4
        hold = 0
        for j = 1:4
            for b = 1:3
                for c = 1:3
                    if epsi(a, b, c) ~= 0 
                        hold = hold + g*epsi(a, b, c)*(A{c}(j))*F_upper_non_abel(j, v, b)
                    end
                end 
            end
        end
        epsi_F_upper_non_abel(v, a) = simplify(hold)
    end    
end 

%% Calculate the current
curr_non_abel = sym('curr_non_abel',[4 3])
for a = 1:3
    for v = 1:4
        hold = 0
        for b = 1:4
            hold = hold + Au{b}(v)*(xid*(T{a}*T{b})*xi + xid*(T{b}*T{a})*xi)
        end
        curr_non_abel(v, a) = ((g^2*rho^2)/8)*simplify(hold)
    end
end

%% Calculate the 4-derivative of non abelian field tensor
co_F_upper_non_abel = sym('co_F_upper_non_abel',[4 3])
for v = 1:4
    for a = 1:3
        co_F_upper_non_abel(v, a) = simplify(D_F_upper_non_abel(v, a) - epsi_F_upper_non_abel(v, a))
    end
end

%% calculate equation of motion
EOM_index = sym('EOM_index',[4 3])
for v = 1:4
    for a = 1:3
        EOM_index(v, a) = simplify(D_F_upper_non_abel(v, a) - epsi_F_upper_non_abel(v, a) - curr_non_abel(v, a))
    end
end