
% Programming Assignment- 02

% Dhanesh Kumar

function main(carrier_freq,velocity)

% Input the carrier frequency and velocity
tic;
f_c=str2double(carrier_freq);
v=str2double(velocity);
f_d= f_c*v/(3e8);   % doppler frequency
w_d=2*pi*f_d;
T_c=1/(4*f_d);      % coherence time
T_s=T_c/10;         % sampling interval
L= 100000;          % No. of Samples in the trace
t=0:T_s:L*T_s;
M=128;              % No. of Multipaths, N = 4*M
n=linspace(1,M,M);
theta=-pi+2*pi*rand(1,1);
alpha=((2*pi.*n)-pi+theta)./(4*M);
psi=-pi+2*pi*rand(1,M);
phi=-pi+2*pi*rand(1,1);
Xc=zeros(1,length(t));
Xs=zeros(1,length(t));
% Construction of trace
for i=1:M
    Xc=Xc+sqrt(2/M)*(cos(psi(i))*cos(w_d*t*cos(alpha(i))+phi));
    Xs=Xs+sqrt(2/M)*(sin(psi(i))*cos(w_d*t*cos(alpha(i))+phi));
end
X=Xc+1i*Xs;
modX=abs(X);
modX2=modX.^2;


% Plot of f(Xc) v/s Xc
figure(1);
[a1,b1]=ksdensity(Xc);
plot(b1,a1,'b',LineWidth=1.5);
hold on;
Xc_th=pdf('Normal',b1,0,1/sqrt(2));
plot(b1,Xc_th,'r',LineWidth=1.5);
histogram(Xc,b1,'Normalization','pdf');
hold off;
legend('Simulated(ks.)','Theoritical','Simulated(hist)');
grid on;
title('PDF of Xc');
xlabel('Xc ---> ');
ylabel('f(Xc) --->');


% Plot of f(|X|) v/s |X|
figure(2);
[a2,b2]=ksdensity(modX,'support','positive');
plot(b2,a2,'b',LineWidth=1.5);
hold on;
modX_th=raylpdf(b2,1/sqrt(2));
plot(b2,modX_th,'r',LineWidth=1.5);
histogram(modX,b2,'Normalization','pdf');
hold off;
legend('Simulated(ks.)','Theoritical','Simulated(hist)');
grid on;
title('PDF of |X|');
xlabel('|X| ---> ');
ylabel('f(|X|) --->');


% Plot of f(|X|^2) v/s |X|^2
figure(3);
[a3,b3]=ksdensity(modX2,'support','positive');
plot(b3,a3,'b',LineWidth=1.5);
hold on;
modX2_th=exp(-b3);
plot(b3,modX2_th,'r',LineWidth=1.5);
histogram(modX2,b3,'Normalization','pdf');
hold off;
legend('Simulated(ks.)','Theoritical','Simulated(hist)');
grid on;
title('PDF of |X|^2');
xlabel('|X|^2 ---> ');
ylabel('f(|X|^2) --->');


% Plot of F(|X|) v/s |X|
figure(4);
[a4,b4]=ecdf(modX);
plot(b4,a4,'b',LineWidth=1.5);
hold on;
modX_cdf_th=raylcdf(b4,1/sqrt(2));
plot(b4,modX_cdf_th,'r',LineWidth=1.5);
hold off;
legend('Simulated','Theoritical');
grid on;
title('CDF of |X|');
xlabel('|X| ---> ');
ylabel('F(|X|) --->');


% Plot of F(|X|^2) v/s |X|^2
figure(5);
[a5,b5]=ecdf(modX2);
plot(b5,a5,'b',LineWidth=1.5);
hold on;
modX2_cdf_th=1-exp(-b5);
plot(b5,modX2_cdf_th,'r',LineWidth=1.5);
hold off;
legend('Simulated','Theoritical');
grid on;
title('CDF of |X|^2');
xlabel('|X|^2 ---> ');
ylabel('F(|X|^2) --->');


% Plot of Auto-correlation of Xs v/s noramlized time(f_d*t)
figure(6);
R_XsXs=autocorr(Xs,length(t)-1);
plot(f_d*t,R_XsXs,'b',LineWidth=1.5);
hold on;
R_XsXs_th=besselj(0,w_d*t');
plot(f_d*t,R_XsXs_th,'r',LineWidth=1.5);
hold off;
legend('Simulated','Theoritical');
axis([0 f_d*t(650) -1.5*max(R_XsXs) 1.5*max(R_XsXs)]);
grid on;
title('Autocorrelation of Xs');
xlabel('Normalized Time: f_d*t ---> ');
ylabel('R_X_s_X_s --->');


% Plot of Auto-correlation of X v/s noramlized time(f_d*t)
figure(7);
R_XX=2*autocorr(X,length(t)-1);
plot(f_d*t,R_XX,'b',LineWidth=1.5);
hold on;
R_XX_th=2*besselj(0,w_d*t');
plot(f_d*t,R_XX_th,'r',LineWidth=1.5);
hold off;
legend('Simulated','Theoritical');
axis([0 f_d*t(650) -1.5*max(R_XX) 1.5*max(R_XX)]);
grid on;
title('Autocorrelation of X');
xlabel('Normalized Time: f_d*t ---> ');
ylabel('R_X_X --->');


% Plot of Cross-correlation of Xc and Xs v/s noramlized time(f_d*t)
figure(8);
R_XcXs=crosscorr(Xc,Xs,length(t)-1);
len=(length(R_XcXs)-1)/2;
plot(f_d*t,R_XcXs(len+1:end),'b',LineWidth=1.5);
hold on;
R_XcXs_th=zeros(length(t),1);
plot(f_d*t,R_XcXs_th,'r',LineWidth=1.5);
hold off;
legend('Simulated','Theoritical');
axis([0 f_d*t(650) -1.5*max(R_XcXs) 1.5*max(R_XcXs)]);
grid on;
title('Crosscorrelation of Xc & Xs');
xlabel('Normalized Time: f_d*t ---> ');
ylabel('R_X_c_X_s --->');

toc;