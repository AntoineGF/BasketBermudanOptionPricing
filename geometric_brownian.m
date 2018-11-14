% geometric_brownian(N,r,sigma,T) simulates a geometric Brownian motion 
% on [0,T] using N normally distributed steps and parameters r and sigma

function [X] = geometric_brownian(N,r,sigma,T)

t = (0:1:N)'/N;                   % t is the column vector [0 1/N 2/N ... 1]

W = [0; cumsum(randn(N,1))]/sqrt(N); % S is running sum of N(0,1/N) variables

t = t*T;
W = W*sqrt(T);

Y = (r-(sigma^2)/2)*t + sigma * W;

X = exp(Y);

% plot(t,X);          % plot the path
% hold on
% plot(t,exp(r*t),':');
% axis([0 T 0 max(1,exp((r-(sigma^2)/2)*T+2*sigma))])
% title([int2str(N) '-step geometric Brownian motion and its mean'])
% xlabel(['r = ' num2str(r) ' and sigma = ' num2str(sigma)])
% hold off

end