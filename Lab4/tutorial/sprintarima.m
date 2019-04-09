% Milan Vojnovic, Jean-Yves Le Boudec, June 04 2003
% Slavisa Sarafijanovic, May 2004.
% Corrected and updated by Alaeddine El Fawal, April 2005
% ARIMA model for Sprint data, NEW

% Acknowledgement:
%  This script was written by taking sarmabat.m as a template, 
%    Hardy Hurd, "Matlab scripts and programs," 
%    http://www.stat.unc.edu/faculty/hurd/stat185Data/progdoc.html
%  The scripts acf.m, acvf.m, and pacf.m originate from: 
%    http://www.math.kth.se/matstat/gru/5b1545/aktuell.htm

clear all
close all
load sprint250.txt
y = sprint250;

nob = length(y)-26; % number of initial data points used for fitting

% samples used for fitting
yfit = y(1:nob);
% the rest is used to compare with forcasts
yval = y(nob+1:length(y)); 

% differenced data at lags 1 and 16
s = 16
xfit = yfit(s+1:length(yfit))-yfit(1:length(yfit)-s);
xfit = diff(xfit);

alpha = 0.05; % for confidence intervals

% we consider ARMA(p+s,q+s), for a p in pp, and q in qq

pp = [1 2 3 4]
qq = [0 1 2 3 4]

aics = zeros(length(pp),length(qq));
baic = 10^6;

for kp=1:length(pp)
 for kq=1:length(qq)
    
    p = pp(kp);
    q = qq(kq);

    disp(sprintf('*** ARMA(%d,%d)\n', p+s, q+s))

    % A is AR polynomial, of the form
    % A(B) = 1 + a_1 B + ... + a_p B^p + a_s B^s + ... + a_{s+p} B^{s+p}
    A=[ones(1,p+1) zeros(1,s-p-1) ones(1,p+1)]; 

    % M is MA polynomial, of the form
    % M(B) = 1 + m_1 B + ... + m_p B^q + m_s B^s + ... + m_{s+q} B^{s+q}
    M=[ones(1,q+1) zeros(1,s-q-1) ones(1,q+1)];

    mod = idpoly(A,[],M,[]);

    % fixedp are indices of non-zero elements of the vector 
    % [A(2) A(3) ... A(length(A)) M(2) M(3) ... M(length(M))]

    zA = find(mod.a==0);
    fixedp = []; 
    if(length(zA))
     fixedp = [fixedp zA-1];
    end
    zM = find(mod.c==0);
    if(length(zM))
     fixedp = [fixedp zM-1+mod.na];
    end

    % set initial values for the coefficients
    A(2:length(A)) = 0.01 * A(2:length(A));
    M(2:length(M)) = 0.01 * M(2:length(M));

    mod = idpoly(A,[],M,[]);

    % fix those coefficients with indices in fixedp to zero
    set(mod, 'FixedParameter', fixedp);

    % estimate coefficients
    %%%% Alaeddine: Error line
    % mod2 = armax(yfit,mod)
    % The correct line is
    mod2 = armax(xfit,mod)

    aki = aic(mod2)
    aics(kp,kq) = aki;

    if(aki < baic) 
     baic = aki;
     baicmodel = mod2;
     baicp = [p q];
    end

 end
end

% Alaeddine: Erroneous line
%disp(sprintf('Best AIC: ARMA(%d,%d)\n', baicp(1)+s, baicp(2)+s))
% The correct Line is:
disp(sprintf('Best AIC: ARMA(%d,%d)\n', baicp(1), baicp(2)))

% compute forecasts of the differenced data
xfore = zeros(length(y)-nob,1);

for k=1:(length(y)-nob)
 xftemp = predict(baicmodel, [xfit; zeros(k,1)], k);
 xfore(k) = xftemp(length(xftemp));
end



% forcasts of the non-differenced data
yf = [yfit; zeros(length(y)-nob,1)];
for k=1:length(y)-nob
    yf(nob+k)=xfore(k)+yf(nob+k-1)+yf(nob+k-s)-yf(nob+k-(s+1));
end

% compute variance of the forecasts, for confidence intervals
syms xx
fa = poly2sym(baicmodel.a(length(baicmodel.a):-1:1), xx);
fm = poly2sym(baicmodel.c(length(baicmodel.c):-1:1), 'xx');
h = length(y)-nob;
%%%%%%% Alaeddine: Error
%csym = taylor(fm/fa/(1-xx)/(1-xx^s), h);
% The correct line is
csym = taylor(fm/(fa*(1-xx)*(1-xx^s)), h);
csym = vpa(csym, 10);
c = sym2poly(csym);

%%%%%%%%%% compute forecasts of the differenced and non-differenced data: method 2
% notation: arma model defined by matrices A and B: A*output = B*whiteNoise
% the data fitting to arma model results in A=baicmodel.a, B=baicmodel.c

A=baicmodel.a;
B=baicmodel.c;

% MA(infinity) representation of the output: output = C*eps
%% in our case: 
% output =  x , xfit, xfore; 
% C is implse response of a filter defined by B and A ; 
% eps is obtained by inverse filtering of x
C=impz(B,A,length(y));
eps=filter(A,B,xfit);

% prediction of the differenced data (see lecture notes)

sigma=std(eps);
xfore2 = zeros(length(y)-nob,1);
mse2 = zeros(length(y)-nob,1);

for k=1:(length(y)-nob)
    xfore2(k)=0;
    for m=1:length(eps)
        xfore2(k) = xfore2(k) + C(k+m)*eps((length(eps)-m+1));
    end
end

mse2(1)=C(1)^2*sigma^2;
for k=2:(length(y)-nob)
    mse2(k)=mse2(k-1)+C(k)^2*sigma^2;
end

% prediction of the non-differenced data

conf2 = norminv(1-alpha/2,0,1)*sqrt(mse2);
yf2 = [yfit; zeros(length(y)-nob,1)];
for k=1:length(y)-nob
    yf2(nob+k)=xfore2(k)+yf2(nob+k-1)+yf2(nob+k-s)-yf2(nob+k-(s+1));
end

% prediction of the variances and confidance intervals
%differencing filter
A3=1;
B3=zeros(1,s+1);
B3(1)=1;
B3(2)=-1;
B3(s)=-1;
B3(s+1)=1;

% MA(infinity) representation: data=C3*eps3
eps3=eps;
C3=filter(A3,B3,C);

sigma3=std(eps3);
mse3 = zeros(length(y)-nob,1);

mse3(1)=C3(1)^2*sigma3^2;
for k=2:(length(y)-nob)
    mse3(k)=mse3(k-1)+C3(k)^2*sigma3^2;
end

conf3 = norminv(1-alpha/2,0,1)*sqrt(mse3);

% plots

figure

plot(200:nob, y(200:nob), 'k')
hold on
plot(nob+1:length(y), y(nob+1:length(y)), 'k')
plot(nob+1:length(yf), yf(nob+1:length(yf)), 'ro-')

sumc2 = cumsum((c(length(c):-1:1).^2)');

conf = norminv(1-alpha/2,0,1)*sqrt(baicmodel.NoiseVariance*sumc2);

plot(nob+1:length(yf), yf(nob+1:length(yf))+conf, 'b')
plot(nob+1:length(yf), yf(nob+1:length(yf))-conf, 'b')
hold off
xlabel(sprintf('Best AIC model: ARMA(p,q) = (%d,%d)', baicp(1)+s, baicp(2)+s));
ylabel('data and forecast')

figure

plot(200:nob, y(200:nob), 'k')
hold on
plot(nob+1:length(y), y(nob+1:length(y)), 'k')
plot(nob+1:length(yf2), yf(nob+1:length(yf2)), 'ro-')

plot(nob+1:length(yf2), yf(nob+1:length(yf2))+conf3, 'b')
plot(nob+1:length(yf2), yf(nob+1:length(yf2))-conf3, 'b')
hold off
xlabel(sprintf('Best AIC model: ARMA(p,q) = (%d,%d)', baicp(1)+s, baicp(2)+s));
ylabel('data and forecast')

%%%%% Alaeddine: Monte Carlo Method to compute the confidence interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samp_nb = 1000;
sigma4 = baicmodel.NoiseVariance;
resids = eps;
lpred = 26; % nb of poins to predict

for i = 1:samp_nb
    resids_fore = normrnd(0,sqrt(sigma4),1,lpred);
    resids_tot = [resids' resids_fore];
    x_inv = filter(B,A,resids_tot);
    x_rlz = x_inv(end - lpred + 1 : end);
    
    xfore_rlz = x_rlz;
    
    y_rlz(i,1:s+2) = yfit(end -s -1 :end)';
    for k=1:lpred
        y_rlz(i,s+2+k)=xfore_rlz(k)+y_rlz(i,s+2+k-1)+y_rlz(i,s+2+k-s)-y_rlz(i,s+2+k-(s+1));
    end
end

y_rlz = y_rlz(:,end-lpred+1:end);

for i = 1:lpred
    T_tmp = sort(y_rlz(:,i));
    cnf_i(i) = T_tmp(26);
    cnf_s(i) = T_tmp(975);
end


figure

plot(200:nob, y(200:nob), 'k')
hold on
plot(nob+1:length(y), y(nob+1:length(y)), 'k')
plot(nob+1:length(yf2), yf(nob+1:length(yf2)), 'ro-')

plot(nob+1:nob+lpred,cnf_i);
plot(nob+1:nob+lpred,cnf_s);

title('Monte-Carlo estimation of the confidence interval');

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












% check residuals

resids = pe(baicmodel,xfit);
racf = xcorr(resids, 'coeff');
racf = racf(ceil(length(racf)/2):length(racf));
racf = racf(1:30);
racvf = acvf(resids);

% compute Ljung-Box statistic and p-values
racf = racvf/racvf(1);
n = length(xfit);

nlags = 10:15;
p=zeros(length(nlags),1);

for k = 1:length(nlags)
  Q = n*(n+2)*sum(racf(2:nlags(k)+1).^2./(n-(1:nlags(k))));
  p(k) = 1-chi2cdf(Q,nlags(k)); 
end

% plot residuals, acf, pacf, p-values of Ljung-Box

figure

subplot(4,1,1)
stem(resids, 'k.-')
ylabel('residuals')

subplot(4,1,2)

acf(resids,1); 
ylabel('ACF of residuals')
axis([0 30 -1 1])

subplot(4,1,3)
pacf(racvf, 30, 1, 1);
ylabel('PACF of residuals')

subplot(4,1,4)
plot(nlags, p, 'k.');
title('p-values of Ljung-Box Chi-Squared Statistics')
ylabel('p-values')



