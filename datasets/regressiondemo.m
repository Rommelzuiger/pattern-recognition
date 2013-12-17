%
%    REGRESSIONDEMO(YBNDS,SIGMA_X,FUNCY)
%
% Show the construction of probabilistic regression models using a
% Bayesian approach (assuming a prior on Y and a generative model
% P(X|Y)).  The prior on Y is a uniform distribution between YBNDS(1)
% and YBNDS(2). The default is YBNDS = [0.5 4.0].
% The generative model P(X|Y) is a Normal distribution where the
% mean is a function of y and the standard deviation is given by
% SIGMA_X. The default is  FUNCY = @(y)sqrt(y) .
% 
% Given the prior and the generative model, the MMSE (minimum mean
% square estimator), the MAP (maximum a posteriori) and the ML (maximum
% likelihood estimator are returned.
function regressiondemo(ybnds,sigma_x,funcy)

if nargin<3
	funcy = @(x)sqrt(x);
end
if nargin<2
	sigma_x = 0.2;
end
if (nargin<1) | isempty(ybnds)
	ybnds = [0.5 4.0];
end

% define the prior on y:
global a b;
a = ybnds(1); b = ybnds(2);
p_y = @uniformprior;

% define the conditional probability p(x|y)
px_y = @(x,y)normalpdf(x,funcy(y),sigma_x);

% define the range of x and y:
h = 0.01;
x = [0.0:h:2.5]';
y = [0.0:h:8]';
V = [0 2.5 0 5];

% find the full distribution:
N = length(y);
res_pxy = zeros(N,length(x));
for i=1:N
	res_px_y(i,:) = px_y(x,y(i));
end
% apply the prior on y:
for i=1:N
	res_pxy(i,:) = px_y(x,y(i))*p_y(y(i));
end
% compute the distribution over x:
px = sum(res_pxy,1)*h;
% finally the posterior:
res_py_x = res_pxy./repmat(px,length(y),1);
% now the solutions for the different estimators:
sol = zeros(3,length(x));
for i=1:length(x)
	sol(1,i) = res_py_x(:,i)'*y(:)*h;
	[maxval,j] = max( res_py_x(:,i) );
	sol(2,i) = y(j);
	[maxval,j] = max( res_px_y(:,i) );
	sol(3,i) = y(j);
end

% now plot various things;
figure(1); clf;
subplot(3,2,1); plot(y,p_y(y));
V0 = axis; V0(1:2)=V(3:4); axis(V0);
xlabel('y'); ylabel('p(y)'); title('Prior of y');

choosen_y = 2.0;
subplot(3,2,2); plot(x,px_y(x,choosen_y));
xlabel('x'); ylabel('p(x|y)');
title(sprintf('p(x|y) for a single y=%5.2f',choosen_y));

subplot(3,2,3); contourf(x,y,res_px_y);
axis(V);
xlabel('x'); ylabel('y'); title('p(x|y) for all y''s');

subplot(3,2,4); contourf(x,y,res_pxy);
axis(V);
xlabel('x'); ylabel('y'); title('p(x,y)=p(x|y)p(y)');

subplot(3,2,5); plot(x,px);
xlabel('x'); ylabel('p(x)'); title('p(x)=\Sigma_y p(x|y)p(y)');

subplot(3,2,6); contourf(x,y,res_py_x,20);
axis(V);
xlabel('x'); ylabel('y'); title('p(y|x)');

figure(2); clf;
contourf(x,y,res_py_x,20);
axis(V);
xlabel('x'); ylabel('y'); title('p(y|x)');
hold on;
h1=plot(x,sol(1,:),'r');
h2=plot(x,sol(2,:),'g');
h3=plot(x,sol(3,:),'y');

legend([h1;h2;h3],'MMSE','MAP','ML','Location','best');

return

function py = uniformprior(y)

global a b;
if isempty(a), a = 0.5; end
if isempty(b), b = 4.0; end

py = (y>=a & y<=b)./(b-a);

return


function px = normaldistr(x,y)
global sigma_x;
if isempty(sigma_x), sigma_x = 0.2; end

if length(y)>1
	error('The y in normaldistr can only have length 1.');
end

px = normalpdf(x,sqrt(y),sigma_x);

return
