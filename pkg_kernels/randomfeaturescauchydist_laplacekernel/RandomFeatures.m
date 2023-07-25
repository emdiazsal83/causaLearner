rand('seed',123)
randn('seed',123)
N = 2e3;
d = 1;
x = randn(N,d);
e = .1*randn(N,d);
y = sin(x) + e;
%
sigma = 1;
% K     = rbf(x,x,sigma);
K = exp(-mandist(x')/sigma);
alpha = (K+.1*eye(N))\y;
yhat  = K*alpha;
%
figure(1),clf,
plot(x,y,'.')
hold on, grid on
plot(x,yhat,'.')
%% Random Fourier Features
D   = 500;
% W   = (1/sigma)*randn(D,d);
m = 0;
bb = 1/sigma;
W = m+bb*tan(pi*(rand(D,d)-1/2));
% mapping cos
b  = 2*pi*rand(D,1);
mc = sqrt(2/D)*cos(x*W' + repmat(b',N,1));
Kc = mc*mc';
% mapping sin cos
msc = (1/sqrt(D))*[sin(x*W') cos(x*W')];
Ksc = msc*msc';
% mapping complex-exponential
mce = (1/sqrt(D))*exp(1i*(x*W'));
Kce = real(mce*mce');
%% orthogonal random features
dD     = max([D,d]);
W      = randn(dD);
G      = W(1:dD,1:dD);
[Q, ~] = qr(G);
sll    = rng;
% rng('default') % seed reproducibility
S      = diag(sqrt(chi2rnd(dD,[1 dD])));
rng(sll);
W      = S*Q/sigma;
W      = W(1:D, 1:d);  % added
% mapping orthogonal feat. complex-exponential
moce   = (1/sqrt(D))*exp(1i*(x*W'));
Koce   = real(moce*moce');
%% plot
figure(2), clf
plot(Kc(:),K(:),'.') % cosinus mapping
hold on, grid on 
plot(Ksc(:),K(:),'.') % sinus cosinus mapping
plot(Kce(:),K(:),'.') % complex exponential mapping
plot(Koce(:),K(:),'.') % orthogonal random features
plot([0 1],[0 1],'-')
legend('cos(wx+b)','[cos, sin]','exp(iwx)','orthog','y=x')
%%
ec   = sqrt(mean((Kc(:)-K(:)).^2));
esc  = sqrt(mean((Ksc(:)-K(:)).^2));
ece  = sqrt(mean((Kce(:)-K(:)).^2));
eoce = sqrt(mean((Koce(:)-K(:)).^2));
figure(3),clf
c = categorical({'cos(wx+b)','[cos, sin]','exp(iwx)','orthog'});
bar(c,[ec esc ece eoce])
set(gca,'yscale','log'), grid on
ylabel('RMSE')

%% Cauchy
% m = 0;
% b = 1/sigma;
% ca = m+b*tan(pi*(rand(D,d)-1/2));
