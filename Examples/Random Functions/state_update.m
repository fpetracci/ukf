function x_new = state_update(x_old,u)
% State update function. 
% Change states update definition to represent your own state function:
% x_dot = f(x) 
% in discrete time.

% states pre allocation
x_new = zeros(size(x_old));

% states update
x_new(1) = x_old(1)^(0.99);
x_new(2) = 1 + sqrt(x_old(1)*x_old(2));

if size(x_old,1) > 2
	for i=3 : size(x_old,1)
		x_new(i) = x_old(i);
	end
end

end

