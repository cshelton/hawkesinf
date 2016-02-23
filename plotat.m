function plotat(pn,b,minb,maxb,nb)

pn = pn/sum(pn);
b = b/sum(b);

delb = (maxb-minb)/nb;
mnbs = minb+(0:(nb-1))*delb;
mxbs = mnbs+delb;

subplot(2,1,1);
plot(pn);
axis([0 14 0 0.5]);
subplot(2,1,2);
plot((mnbs+mxbs)/2,b);
