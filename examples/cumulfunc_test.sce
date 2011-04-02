format(20,20);
x=0.6;
y=0.4;
a=2.0;
b=3.0;
[p,q] = cdfbet("PQ",x,y,a,b);
printf ("cdfbet : p = %.18f\tq = %.18f\n", p, q);

s=2;
xn=5;
pr=0.7;
ompr=0.3;
[p,q] = cdfbin("PQ",s,xn,pr,ompr);
printf ("cdfbin : p = %.18f\tq = %.18f\n", p, q);

p=0.7;
q=0.3;
x=2;
df = cdfchi("Df", p,q,x);
printf ("cdfchi : df = %.18f\n", df);
df = 6;
x = cdfchi("X",df,p,q);
printf ("cdfchi : x = %.18f\n", x);

p=0.1;
q=0.9;
x=2;
df=2;
pnonc = cdfchn("Pnonc",p,q,x,df);
printf ("cdfchn : pnonc = %.18f\n", pnonc);

p=0.1;
q=0.9;
f=2;
dfn=3;
dfd = cdff("Dfd",p,q,f,dfn);
printf ("cdff : dfd = %.18f\n", dfd);


p=0.1;
q=0.9;
f=2;
dfn=3;
dfd=5;
pnonc = cdffnc("Pnonc",p,q,f,dfn,dfd);
printf ("cdffnc : pnonc = %.18f\n", pnonc);

p=0.1;
q=0.9;
x=2;
shape=3;
scale = cdfgam("Rate",p,q,x,shape);
printf ("cdfgam : scale = %.18f\n", scale);

s=2.0;
xn=3.0;
pr=0.7;
ompr=0.3;
[p,q] = cdfnbn("PQ",s,xn,pr,ompr);
printf ("cdfnbn : p = %.18f\tq = %.18f\n", p, q);

x=2.0;
mean=0.0;
sd=1.0;
[p,q] = cdfnor("PQ",x,mean,sd);
printf ("cdfnor : p = %.18f\tq = %.18f\n", p, q);

p=0.4;
q=0.6;
s=5.0;
xlam = cdfpoi("Xlam",p,q,s);
printf ("cdfpoi : xlam = %.18f\n", xlam);

p=0.4;
q=0.6;
t=-5.0;
df = cdft("Df",p,q,t);
printf ("cdft : df = %.18f\n", df);

