# ukf

The filter needs two external function that the user has to write:

state_update(x,u)	- how the state is updated given state and input
state2meas(x,u)		- how the measurement vector is obtained given state and input

