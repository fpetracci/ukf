function y = state2meas(x, u)
% state2meas function that calculates outputs given state.

	global C_td D_td
	
    y = C_td * x + D_td * u ;
end