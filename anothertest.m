function [pn,b] = anothertest(p0,pi,mu,sigma,m,minb,maxb,nb)

	x = [0];
	y = randn(1)*sigma+mu;
	pn = [];

	b = zeros(nb,1);
	delb = (maxb-minb)/nb;
	mnbs = minb+(0:(nb-1))*delb;
	mxbs = mnbs+delb;

	mu0 = -2;

	q = [0.3 0.2 0.5];

	for i=1:m
		if (i>m/2)
			if (length(x)>length(pn))
				pn(length(x)) = 1;
			else
				pn(length(x)) = pn(length(x))+1;
			end
			ys = sum(y);
			b(mnbs<=ys & mxbs>ys) = b(mnbs<=ys & mxbs>ys) +1;
		end;
		rr = rand(1);
		if (rr<q(1))
			j = ceil(rand(1)*length(x));
			y(j) = randn(1)*sigma+mu;
		else
			rr = rr-q(1);
			if (rr<q(2))
				if (length(x)==1)
					ratio = q(3)*p0*(1-pi)/(q(2)*(1-p0));
				else
					ratio = q(3)*pi/q(2);
				end;
				ratio = ratio*(1/(sqrt(2*pi)*sigma)*exp(-(mu-mu0)*(mu-mu)/sigma/sigma/2));
				if (rand(1)<ratio) 
					x(length(x)) = 1;
					x(length(x)+1) = 0;
					y(length(x)) = randn(1)*sigma+mu;
					%y(length(x)) = mu0;
				end;
			else
				if (length(x)==2)
					ratio = q(3)*p0*(1-pi)/(q(2)*(1-p0));
				else
					ratio = q(3)*pi/q(2);
				end;
				ratio = ratio*(1/(sqrt(2*pi)*sigma)*exp(-(mu-mu0)*(mu-mu)/sigma/sigma/2));
				ratio = 1/ratio;
				if (length(x)>1 && rand(1)<ratio) 
					x = x(1:end-1);
					y = y(1:end-1);
				end;
			end;
		end;
	end;
