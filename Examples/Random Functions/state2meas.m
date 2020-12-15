function meas = state2meas(x)
% given state vector, it outputs its virtual measurements vector.
% Change outputs function definition to represent your own function:
% y = h(x) 
% in discrete time.

% outputs vector pre-allocation
n_out = 3;
meas = zeros(n_out, 1);

% outputs function definition
if n_out <= size(x,1)
	for i = 1:n_out
		meas(i) = log(x(i)^2);
	end
else
	for i = 1:size(x,1)
		meas(i) = log(x(i)^2);
	end
	for j = i+1:n_out
		meas(j) = x(j-i);
	end
end
