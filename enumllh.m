function logw = enumllh(ts)

% samples from K(i,j,s,t) = M(i,j)e^(-t-s)
%              M = [1 2; 1 1]/4
%              mu_1 = 1/10
%              mu_2 = 1/5
%
% evidence is known everywhere except for event type 1 from 1 to 3
% no events except at event type 2 at time 4
% time is from 0 to 5

M = [1 2; 1 1]/4;
t0 = 1;
t1 = 3;
T = 5;
etime = 4;
mu = [0.1 0.0001];

logw = -T*mu(1)-T*mu(2);
logw = logw-Kbar(M,2,2,etime,etime,T)-Kbar(M,2,1,etime,etime,T);
ew = mu(2);
for t=ts,
	logw = logw-Kbar(M,1,1,t,t,T)-Kbar(M,1,2,t,t,T);
	ew = ew+K(M,1,2,t,etime);
	
	g = mu(1);
	for s=ts
		if (s<t)
			g = g+K(M,1,1,s,t);
		else
			break;
		end;
	end;
	logw = logw+log(g);
end;
logw = logw+log(ew);

function k = K(M,i,j,s,t)
	if (s>t)
		k = 0;
	else
		k = M(i,j)*exp(-(t-s));
		%k = M(i,j)/((t-s+1)*(t-s+1));
	end;

function p = Kbar(M,i,j,s0,s,t)
	if (s>t)
		p = 0;
	else
		p = M(i,j)*(exp(s0-s)-exp(s0-t));
		%p = M(i,j)*(1/(s-s0+1) - 1/(t-s0+1));
	end;
