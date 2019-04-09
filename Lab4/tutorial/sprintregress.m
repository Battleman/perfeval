% Milan Vojnovic, May 21, 2003
% Slavisa Sarafijanovic, May 24, 2004
% Alaeddine El Fawal, April 2005
% Corrected by Slavisa Sarafijanovic, May 29, 2006

close all;
clear all;

load sprint250.txt
y = sprint250;

alpha = 0.95;

% we assume we do not know last 26 data points
nob = 224
t = [1:nob]';

X=zeros(nob, 7);

X(:,1)=1;
X(:,2)=t;
X(:,3)=t.^2;
X(:,4)=cos(2*pi*t/16);
X(:,5)=sin(2*pi*t/16);
X(:,6)=cos(2*pi*t/8);
X(:,7)=sin(2*pi*t/8);

[b,bint,r,rint,stats] = regress(y(1:nob), X, alpha);

t= [1:length(y)]';

yreg=b(1)+b(2)*t+b(3)*t.^2+b(4)*cos(2*pi*t/16)+b(5)*sin(2*pi*t/16);
yreg=yreg+b(6)*cos(2*pi*t/8)+b(7)*sin(2*pi*t/8);

figure
plot(t(1:nob), y(1:nob), 'k')
hold
plot(t(1:nob), yreg(1:nob), 'r')
hold off
legend('data', 'regression')
npast=16;

%err = y(nob+1:length(y))-yreg(nob+1:length(yreg));
res = y(1:nob)-yreg(1:nob);

%cint = tinv((1+alpha)/2, length(err)-1)*std(err)/sqrt(length(err));
cint = norminv((1+alpha)/2)*std(res,1)*(nob-1)/(nob-length(b));
% cint = norminv((1+alpha)/2)*sqrt(sum(res.^2)/(nob-length(b))); <- gives very similar as previous line



%%%%%%%%%%%% added by Slavisa, May 27 2006 (with comments)
%%% computed cint doesn't take correctly into account uncertanty caused by
%%% estimating b(i), i=1,...,7.
%%% we can compute ci - the ci for predicted yreg using (Thm 8.1.1)

% first compute s2 and K using model X of dimension nob*7 for which
% observations exist, and from which b(i) are estimated
ssr=sum((y(1:nob)-yreg(1:nob)).^2);
s2=ssr/(nob-length(b));
K=(X'*X)\X';

% then (periodicity) assume validity of the model extends to nob+1,...,250.
t = [1:length(y)]';
Xe=zeros(length(y), 7);
Xe(:,1)=1;
Xe(:,2)=t;
Xe(:,3)=t.^2;
Xe(:,4)=cos(2*pi*t/16);
Xe(:,5)=sin(2*pi*t/16);
Xe(:,6)=cos(2*pi*t/8);
Xe(:,7)=sin(2*pi*t/8);

% and compute vector of variance bias terms g and conf. int. ci using (Thm 8.1.1.4.c)
g=sum((Xe*K).^2,2);
eta=tinv((1+alpha)/2, nob-length(b));
ci=eta*((g*s2).^0.5);

% we are interested in "prediction with confidence intervals" (lecture notes, page 203)
% which is probably meant to be prediction interval defined as: 
% new sample will be in the prediction interval with probability = alpha 
% it can be computed as in Question 8.1.2 (lecture notes, page 179)
pi=eta*(((1+g)*s2).^(0.5));

% add ci and pi to the figure to compare against cint
figure
plot(t(nob-npast:nob), y(nob-npast:nob), 'k')
hold on
plot(t(nob+1:length(y)), y(nob+1:length(y)), 'o-b')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg)), 'r')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))+cint, 'r-.')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))-cint, 'r-.')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))+ci(nob+1:length(yreg)), 'g-')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))-ci(nob+1:length(yreg)), 'g-')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))+pi(nob+1:length(yreg)), 'g-.')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))-pi(nob+1:length(yreg)), 'g-.')
hold off

%%% conclusion on pi vs cint: in this example, pi and cint are very similar beacuse
%%% the variance bias terms are small; in other words, the uncertanty due to estimating b(i) 
%%% is much smaller then the uncertanty due to the (estimated) variance of residuals
%%% BUT, in general, this might not be true and should be checked.
%%% or simply always use pi (as computed above).

%%%%%%%%%%%%%%%% end of added by slavisa

%plot(t(nob-npast:nob), y(nob-npast:nob), 'o-k')
figure
plot(t(nob-npast:nob), y(nob-npast:nob), 'k')

hold on
plot(t(nob+1:length(y)), y(nob+1:length(y)), 'o-b')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg)), 'r')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))+cint, 'r-.')
plot(t(nob+1:length(yreg)), yreg(nob+1:length(yreg))-cint, 'r-.')
hold off


