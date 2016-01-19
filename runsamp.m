function [ans,wts] = runsamp(n)

s = [0,0];
swt = 0;
lws = [];
ahist = [];
for i=1:n
	[lw,ts] = simpsamp;
	lws = [lws lw];
	w = exp(lw);
	s(1) = s(1)+sum(ts<2)*w;
	s(2) = s(2)+sum(ts>=2)*w;
	swt = swt+w;
	ans = s/swt;
	ahist = [ahist; ans];
	if (mod(i,10000)==0)
		subplot(1,2,1);
		hist(lws,20);
		subplot(1,2,2);
		plot(ahist);
		drawnow;
	end;
end;
	
