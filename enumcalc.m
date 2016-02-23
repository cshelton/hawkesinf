function ans = enumcalc(dt,maxn)

ttl = 0;
den = 0;

mint = 1+dt/2;
maxt = 3+dt/2;

vol = 1;
for n=0:maxn,
	ts = zeros(1,n);
	for i=1:n
		ts(i) = mint+(i-1)*dt;
	end;
	while (1),
		wt = exp(enumllh(ts))*vol;
		%[wt ts]
		den = den + wt;
		ttl = ttl + wt*n;
		for i=n:-1:0
			if (i==0)
				break;
			end;
			ts(i) = ts(i)+dt;
			if (i==n && ts(i)<maxt)
				break;
			end;
			if (i<n && ts(i)<ts(i+1))
				break;
			end;
			if (i==1)
				ts(i) = mint;
			else
				ts(i) = ts(i-1)+dt*2;
				if (i==n && ts(i)>=maxt)
					ts(i)=maxt;
				end;
				if (i<n && ts(i)>=ts(i+1))
					ts(i) = ts(i+1)-dt;
				end;
			end;
		end;
		if (i==0)
			break;
		end;
	end;
	vol = vol*dt;
	[n den ttl/den]
end;
ttl
den
ans = ttl/den;
