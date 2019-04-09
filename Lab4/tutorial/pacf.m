function a=pacf(g,h,plott,cl)
% pa= pacf(g) gives the Partial ACF for the ACVF g
% pa=pacf(g,h) gives Partial ACF for a ACVF g up to lag h
% pa=pacf(g,h,plott) where plott arbitrary number also gives a plott
% pa=pacf(g,h,plott,cl) draws bounds +- 1.96/sqrt(n)
% if second argument is , the h=length(g),
% note that a second argument must be given, if the plott shall be drawn
% Brockwell page 69 and 92

if nargin==1 | isempty(h), h=length(g); end
f=g(2)/g(1);
v=g(1)*(1-f^2);
a(1)=1;
a(2)=f;
for k=3:h
        fn=(g(k)-f'*g(k-1:-1:2)')/v;
        f=f-fn*f(length(f):-1:1);
        f=[f;fn];
        v=v*(1-fn^2);
        a(k)=fn;
 end
 
 if nargin>=3
    stem(a, 'k.-')
 end
 
 if nargin==4
    hold on
    plot([0 h],1.96*[1 1]/sqrt(length(g)),':')
    plot([0 h],-1.96*[1 1]/sqrt(length(g)),':')
    hold off
 end
 
