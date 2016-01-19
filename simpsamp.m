function [logw,ts] = simpsamp

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

ts = [];
% sample base events:
t = t0;
lastts = [];
while(t<t1)
	t = t-log(rand(1,1))/mu(1);
	if (t<t1)
		lastts = [lastts t];
	end;
end;
while(length(lastts)>0)
	nextts = [];
	ts = [ts lastts];
	for tt=lastts
		torig = tt;
		t = tt;
		while(t<t1)
			t = torig+Kbarinv(M,1,1,t-torig,-log(rand(1,1)));
			if (t<t1)
				nextts = [nextts t];
			end;
		end;
	end;
	lastts = nextts;
end;

%ts = [1.63708 2.94468];

logw = -(T-(t1-t0))*mu(1)-T*mu(2);
ew = mu(2);
for t=ts,
	logw = logw-Kbar(M,1,1,t,t1,T)-Kbar(M,1,2,t,t,T);
	ew = ew+K(M,1,2,t,etime);
end;
logw = logw+log(ew);

function k = K(M,i,j,s,t)
	if (s>t)
		k = 0
	else
		k = M(i,j)*exp(-(t-s));
	end;

function p = Kbar(M,i,j,s0,s,t)
	if (s>t)
		p = 0;
	else
		p = M(i,j)*(exp(s0-s)-exp(s0-t));
	end;

function t = Kbarinv(M,i,j,tcurr,r)
	tr = r/M(i,j)*exp(tcurr);
	if (tr>=1)
		t = Inf;
	else
		t = tcurr-log1p(-tr);
	end;
