function y=acf(x,plott)
% y=acf(x,plot) gives autocorrelation function
% if second argument is given, plot is drawn and appr 95% limits

y=acvf(x);
y=y/y(1);
if nargin ==2
   n=length(y);
   z=zeros(1,n);
   x=[0:(n-1);0:(n-1)];
   a=[z;y];
   plot(x,a,'k');
   hold on;
   plot([0 n-1],[1.96/sqrt(n) 1.96/sqrt(n)],':b',[0 n-1],-ones(1,2)*1.96/sqrt(n),':b');
   hold off;
end

