function y = form(N,l,itt)

y = zeros(size(l));

for n=0:itt
	y = y + N*exp(-l).*(l.^n)/(N+n)/factorial(n);
end;
