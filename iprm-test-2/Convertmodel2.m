function [ c, A, b, free_variables, objective_const_term ] = Convertmodel2( c, A, b, lb, ub, sense )

n = size(A,2);
m = size(A,1);
if (size(lb,2) > 1)
    lb = lb';
end
if (size(ub,2) > 1)
    ub = ub';
end
if (size(sense,2) > 1)
    sense = sense';
end
if (size(c,2) > 1)
    c = c';
end
if (~issparse(A))
    A = sparse(A);
end

if (size(c,1) ~= n || size(sense,1) ~= m || size(lb,1) ~= n || size(ub,1) ~= n)
    error('Incorrect input dimensions')
end

num_of_slacks = 0; 
free_variables = [];
objective_const_term = 0; 
extra_constraints = 0; 
[rows,cols,v] = find(A);
b_new = [];

for i = 1:m    
    if ( sense(i) == '<')          
       num_of_slacks = num_of_slacks + 1;
       rows = [rows; i];
       cols = [cols; n + num_of_slacks];
       v = [v; 1];       
    elseif ( sense(i) == '>')             
       num_of_slacks = num_of_slacks + 1;
       rows = [rows; i];
       cols = [cols;  n + num_of_slacks];
       v = [v; -1]; 
    end   
end

for i = 1:n    
    if ((ub(i) == Inf) && (lb(i)> -Inf))     
        
        if (lb(i) ~= 0)
            b(:) = b(:) - A(:,i).*lb(i);  
            objective_const_term = objective_const_term + c(i)*lb(i);
        end
    elseif ((lb(i) == -Inf) && (ub(i) < Inf))    
        
        k_max = size(cols,1);
        for k = 1:k_max
            if (cols(k) == i)
                v(k) = -v(k);
            end
        end  
        objective_const_term = objective_const_term + c(i)*ub(i);
        c(i) = -c(i); 
        if (ub(i) ~= 0)
            b(:) = b(:) - A(:,i).*ub(i); 
        end
    elseif ((lb(i) > -Inf) && (ub(i) < Inf)) 
       
        if (lb(i) ~= 0)
            b(:) = b(:) - A(:,i).*lb(i);
            objective_const_term = objective_const_term + c(i)*lb(i);
        end
        
        extra_constraints = extra_constraints + 1; 
        num_of_slacks = num_of_slacks + 1;
        b_new = [b_new; ub(i) - lb(i)]; 
        rows = [rows; m + extra_constraints ; m + extra_constraints];
        cols = [cols; i ; n + num_of_slacks];
        v = [v; 1; 1]; 
    else          
        free_variables = [free_variables; i]; 
    end

end

b = [b; b_new];
c = [c; zeros(num_of_slacks,1)];
A = sparse(rows,cols,v,m+extra_constraints,n+num_of_slacks);
size(free_variables)
A = [A,-A(:,free_variables)];
c = [c;-c(free_variables)];
free_variables = [];
end

