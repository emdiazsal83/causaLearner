% Testing Mooij's hsic

addpath(genpath('/home/soulivanh/Documents/proyectos/indepReg/Mooij/code_mooijs/fasthsic'))


bandwidthX = 1;
bandwidthY = 1;

n = 100;
px = 3;
py = 2;

x = rand(n, px); % continuous data A
%x = unidrnd(10, 1, n);% discrete data A

y = randn(n, py);
%y = unidrnd(10,p,n);

% got to make the c++ version first otherwise nrperm has to be zero coz
% matlab version only has that implemented
nrperm = 0;

%nrperm = |nrperm| is the number of permutations used for estimating the p-value
%                    (if > 0,  use original biased HSIC estimator,
%                     if == 0, use gamma approximation,
%                     if < 0,  use unbiased HSIC estimator)

fasthsic(x, y, bandwidthX, bandwidthY, nrperm)

[p, res] = fasthsic(x, y);

res