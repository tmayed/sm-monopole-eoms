%%%% vect4.m

%% spherical 4 vector calculus functions
function output = vect4(func, input, opt1, opt2)
    % check if optional variables used
    if ~exist('opt1','var')
        opt1 = 0
    end
    if ~exist('opt2','var')
        opt2 = 0
    end
    %controller
    if strcmp(func,"grad")
        output = grad(input, opt1)
    elseif strcmp(func,"grad_i")
        output = grad_i(input, opt1, opt2)
    elseif strcmp(func,"grad_comp")
        output = grad_comp(input, opt1)    
    elseif strcmp(func,"grad_comp_i")
        output = grad_comp_i(input, opt1, opt2)     
    elseif strcmp(func,"div")
        output = div(input, opt1)
    elseif strcmp(func,"div_i")
        output = div_i(input, opt1, opt2)    
    elseif strcmp(func,"lap")
        output = lap(input)
    elseif strcmp(func,"lap4")
        output = lap(input)
    else
        %%something has gone wrong!
        output = 666
    end
end 

%% calculate the gradient 
function output = grad(input, opt1)
    %option 1 : apply metric
    output = [grad_i(input, 1, opt1);grad_i(input, 2, opt1);grad_i(input, 3, opt1);grad_i(input, 4, opt1)]
end

%% gradient function by index (for loops)
function output = grad_i(input, opt1, opt2)
    %option 1 : grad index
    %option 2 : apply metric
    syms t r theta phi
    t_term = grad_comp_i(input, 1, opt2)*[1;0;0;0]
    r_term = grad_comp_i(input, 2, opt2)*[0;sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)]
    th_term = grad_comp_i(input, 3, opt2)*[0;cos(theta)*cos(phi);cos(theta)*sin(phi);-sin(theta)]
    ph_term = grad_comp_i(input, 4, opt2)*[0;-sin(phi);cos(phi);0]
    %t term
    if opt1 == 1
        output = simplify(t_term(1, 1) + r_term(1, 1) + th_term(1, 1) + ph_term(1, 1))
    %r term    
    elseif opt1 == 2
        output = simplify(t_term(2, 1) + r_term(2, 1) + th_term(2, 1) + ph_term(2, 1))
    %theta term
    elseif opt1 == 3
        output = simplify(t_term(3, 1) + r_term(3, 1) + th_term(3, 1) + ph_term(3, 1))
    %phi term    
    elseif opt1 == 4
        output = simplify(t_term(4, 1) + r_term(4, 1) + th_term(4, 1) + ph_term(4, 1))
    end
end

%% calculate gradient componets
function output = grad_comp(input, opt1)
    %option 1 : apply metric
    t_term = grad_comp_i(input, 1, opt1)
    r_term = grad_comp_i(input, 2, opt1)
    th_term = grad_comp_i(input, 3, opt1)
    ph_term = grad_comp_i(input, 4, opt1)
    output = [t_term;r_term;th_term;ph_term]
end

%% calculate individual elements of components of gradient 
function output = grad_comp_i(input, opt1, opt2)
    %option 1 : grad index
    %option 2 : apply metric
    syms t r theta phi
    %t term
    if opt1 == 1
        %check if apply metric
        if opt2 == 1
            output = -simplify(diff(input, t))
        else 
            output = simplify(diff(input, t))
        end
    %r term
    elseif opt1 == 2
        output = simplify(diff(input, r))
    %theta term
    elseif opt1 == 3
        output = simplify((1/r)*diff(input, theta))
    %phi term
    elseif opt1 == 4
        output = simplify((1/(r*sin(theta)))*diff(input, phi))
    end
end

%% calculate divergence of 4-vector
function output = div(input, opt1)
    %option 1 : apply metric
    t_term = div_i(input(1, 1), 1, opt1)
    r_term = div_i(input(2, 1), 2, opt1)
    th_term = div_i(input(3, 1), 3, opt1)
    ph_term = div_i(input(4, 1), 4, opt1)
    output = simplify(t_term + r_term + th_term + ph_term)
end

%% calculate individual terms of divergence 
function output = div_i(input, opt1, opt2)
    %option 1 : divergence terms
    syms t r theta phi
    %t term
    if opt1 == 1
        %check if apply metric
        if opt2 == 1
            output = -simplify(diff(input, t))
        else 
            output = simplify(diff(input, t))
        end
    %r term
    elseif opt1 == 2
        output = simplify((1/r^2)*diff((r^2)*input, r))
    %theta term
    elseif opt1 == 3
        output = simplify((1/(r*sin(theta)))*diff(sin(theta)*input, theta))
    %phi term
    elseif opt1 == 4
        output = simplify((1/(r*sin(theta)))*diff(input, phi))
    end
end

%% calculate the laplacian
function output = lap(input)
    syms t r theta phi
    t_term = -diff(diff(input, t), t)
    r_term = diff((r^2)*diff(input, r), r)/(r^2)
    th_term = diff(sin(theta)*diff(input, theta), theta)/((r^2)*sin(theta))
    ph_term = diff(diff(input, phi), phi)/(r*sin(theta))^2
    output = simplify(simplify(t_term)+simplify(r_term)+simplify(th_term)+simplify(ph_term))
end

%% calculate the laplacian of a 4 vector
function output = lap4(input)
    output = [lap(input(1, 1));lap(input(2, 1));lap(input(3, 1));lap(input(4, 1))]
end