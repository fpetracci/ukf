function  [r,Lf_full,Lg_full] = relative_degree(f,g,y,x)
% funzione che ritorna il grado relativo di un sistema che ha come uscita
% y.
% In uscita ritorna rispettivamente, il grado relativo, e tutte le derivate
% direzionale della funzione scalare y(x) lungo il campo vettoriale f e
% g,DA CORREGERE QUESTA DEFINIZIONE, LEGGERA IMPRECISIONE
Lf_y = jacobian(y,x)*f;
Lg_y = jacobian(y,x)*g;
Lf_full = [Lf_y];
Lg_full = [Lg_y];
r = 1; 
while r<length(x) && Lg_y == 0
	r=r+1;
	Lg_y = jacobian(Lf_y,x)*g;
	Lf_y = jacobian(Lf_y,x)*f;
	Lf_full = [Lf_full Lf_y];
	Lg_full = [Lg_full Lg_y];
end

if r==length(x) && Lg_y == 0
	r = -1;
end
