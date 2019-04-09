% Milan Vojnovic, May 21, 2003

load sprint250.txt
y = sprint250;

alpha = 0.05;

nob = 224;

yob = y(1:nob);

d16yob = yob(16+1:length(yob))-yob(1:length(yob)-16);
d116yob = diff(d16yob);


figure

subplot(2,1,1)

acf(d116yob,1); 

axis([0 100 -0.5 1])
xlabel('lags')
ylabel('ACF')

subplot(2,1,2)

pacf(acvf(d116yob), 100, 1, 1);

xlabel('lags')
ylabel('PACF')

