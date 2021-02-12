function  [r_mimo,Lf_full_mimo, T, E] = relative_degree_mimo(f,G,H,x)
% input
	%f = drift
	%G = matrix of g vectors, each one related to a single input
	%H = output functions h_i(x) forall y_i
	%x = states
	
% output
	%r_mimo = vector of relative degrees forall y_i
	%Lf_full_mimo = Lf_full forall y_i
	
n_out = size(H, 1); %output in riga
n_in = size(G, 2);
r_mimo = zeros(n_out, n_in);
for i=1:n_out
	r_min = length(x);
	clear E_tmp
	for j=1:n_in
		[r_mimo(i, j),Lf_full,Lg_full] = relative_degree(f,G(:,j),H(i,:),x);
		if j==1
			E_tmp(:, j) = transpose(Lg_full);
		end
		if r_mimo(i, j) < r_min && r_mimo(i, j)>=0
			in_r = j;
			Lf_full_min = Lf_full;
			r_min = r_mimo(i, j);
			E_tmp = E_tmp(1:r_min, :);
		end
		E_tmp(:, j) = transpose(Lg_full(1:r_min));
	end
	E(i,:) = E_tmp(r_min, :);
	
	if i==1
		Lf_full_mimo = [H(i,:); transpose(Lf_full_min(:, 1:r_min-1))];
		T = [Lf_full_min(:, r_min)];
	else
		Lf_full_mimo = [Lf_full_mimo; H(i,:); transpose(Lf_full_min(:, 1:r_min-1))];
		T = [T; Lf_full_min(:, r_min)];
	end
end

