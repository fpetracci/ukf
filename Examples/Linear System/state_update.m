function x_new = state_update(x_old, u)
% state_update Local function that calculates how the state updates.

	global A_td B_td

    x_new = A_td * x_old + B_td * u;
end

