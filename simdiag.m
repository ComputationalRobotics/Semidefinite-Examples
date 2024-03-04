function [Q,varargout] = simdiag(varargin)
%
% [Q,D1,...,Dm] = simdiag(A1,...,Am,options);
%	Simultaneously diagonalize all input matrices
%
%	Input:		A1,...,A2	pairwise commuting complex normal matrices
%				options		options field:
%								tol: tolerance
%
%	Ouput:		Q			unitary matrix containing the simultaneous
%							eigenvectors of A1,...,Am
%				D1,...,Dm	eigenvalues of A1,...,Am, respectively
%
%	Copyright (c) 2008-2009, Christian Mendl
%	All rights reserved.

if isstruct(varargin{end})
	tol = varargin{end}.tol;
	varargin = varargin(1:end-1);
else
	tol = eps^1.5;
end

n = size(varargin{1},1);

% preprocessing step
Q = dodo(varargin{:});
for j=1:length(varargin)
	varargin{j} = Q'*varargin{j}*Q;
end

% Reference:
%	Angelika Bunse-Gerstnert, Ralph Byers, and Volker Mehrmann,
%	Numerical Methods for Simultaneous Diagonalization,
%	SIAM J. Matrix Anal. Appl. Vol. 14, No. 4, pp. 927-949, October 1993
calc_off2 = @(A)norm(A-diag(diag(A)),'fro')^2;
max_iter = 1e5;
num_iter = 0;
off2 = 0; nscale = 0;
for m=1:length(varargin)
	off2 = off2 + calc_off2(varargin{m});
	nscale = nscale + norm(varargin{m},'fro');
end
while off2 > tol*nscale
	for j=1:n
		for k=j+1:n
			v = zeros(length(varargin),3);
			for m=1:length(varargin)
				v(m,:) = [varargin{m}(j,j)-varargin{m}(k,k),varargin{m}(j,k),varargin{m}(k,j)];
			end
			[c,s] = approx_min(v);
			Q = timesR(Q,j,k,c,s);
			for m=1:length(varargin)
				varargin{m} = rotate(varargin{m},j,k,c,s);
			end
		end
	end
	off2 = 0; nscale = 0;
	for m=1:length(varargin)
		off2 = off2 + calc_off2(varargin{m});
		nscale = nscale + norm(varargin{m},'fro');
	end
	% number of iterations
	num_iter = num_iter + 1;
	if num_iter > max_iter
		fprintf('Exiting: Maximum number of iterations exceeded. Current relative error: %g.\n',off2/nscale);
		break;
	end
end

% eigenvalues
for j=1:min(length(varargin),nargout-1)
	varargout{j} = varargin{j}.*eye(size(Q));
end

	
%% A = A*R, R = R(j,k,c,s)
function A = timesR(A,j,k,c,s)

% A = A*R
A(:,[j,k]) = [c*A(:,j)+s*A(:,k),-conj(s)*A(:,j)+conj(c)*A(:,k)];


%% A = R'*A*R, R = R(j,k,c,s)
function A = rotate(A,j,k,c,s)

A(:,[j,k]) = [c*A(:,j)+s*A(:,k),-conj(s)*A(:,j)+conj(c)*A(:,k)];	% A = A*R
A([j,k],:) = [conj(c)*A(j,:)+conj(s)*A(k,:);-s*A(j,:)+c*A(k,:)];	% A = R'*A


%%
% Approximate minimizer of
% |s c conj(v(:,1)) - c^2 conj(v(:,2)) + s^2 conj(v(:,3))|^2 + |s c v(:,1) + s^2 v(:,2) - c^2 v(:,3)|^2
% for c = cos(theta), s = exp(i phi) sin(theta)
function [c,s] = approx_min(v)

target = @(c,s,v) norm(s*c*conj(v(:,1))-c^2*conj(v(:,2))+s^2*conj(v(:,3)),2)^2+norm(s*c*v(:,1)+s^2*v(:,2)-c^2*v(:,3),2)^2;

[c,s] = calc_min(v(1,1),v(1,2),v(1,3));
m = target(c,s,v);
for j=2:size(v,1)
	[c1,s1] = calc_min(v(j,1),v(j,2),v(j,3));
	x = target(c1,s1,v);
	if x < m
		m = x;
		c = c1; s = s1;
	end
end


%%
% Exact minimizer of
% |s c conj(a0) - c^2 conj(a21) + s^2 conj(a12)|^2 + |s c a0 + s^2 a21 - c^2 a12|^2;
% Refer to
%	H. H. Goldstine and L. P. Horwitz, A Procedure for the
%	Diagonalization of Normal Matrices, J. ACM (1959)
function [c,s] = calc_min(a0,a21,a12)

u = real(a0);
v = imag(a0);

tmp = (a21+conj(a12))/2;
r = abs(tmp); beta = angle(tmp);

tmp = (a21-conj(a12))/2;
s = abs(tmp); gamma = angle(tmp);

nu = beta-gamma;
sin_nu = sin(nu);
cos_nu = cos(nu);

L = u*v-4*r*s*sin_nu;
M = u^2-v^2+4*(r^2-s^2);

A = L*(r^2-s^2)*sin_nu+M*r*s;
B = L*(r^2+s^2)*cos_nu;
C = L*(r^2-s^2)+M*r*s*sin_nu;

tmp = r*s*cos_nu*sqrt(M^2+4*L^2);
phi = (atan2(-A*C+B*tmp,B*C+A*tmp)-beta-gamma)/2;

r_cos_ba = r*cos(beta+phi);
s_sin_ga = s*sin(gamma+phi);
kappa = u^2+v^2-4*(r_cos_ba^2+s_sin_ga^2);
lambda = 4*(u*r_cos_ba+v*s_sin_ga);
theta = -atan2(-lambda,kappa)/4;

c = cos(theta);
s = exp(1i*phi)*sin(theta);


%% "DODO" (diagonalize one then diagonalize the other) preprocessing step
% note: this is instable in general!
function Q = dodo(varargin)

Q = eye(size(varargin{1}));
IE = [ones(1,length(Q)-1),0];

for m=1:nargin
	A = varargin{m};
	% handle degenerate eigenspaces
	A = Q'*A*Q;
	epsilon = 10*max(size(A))*eps(normest(A));
	for j=find((IE==1).*[1,IE(1:end-1)==0])	% unaffected by 'IE' update
		k = j+find(IE(j:end)==0,1,'first')-1;
		[V,d] = schur(A(j:k,j:k)); d = diag(d);
		% handle case d = [1, i, (1-eps)*i]
		[t,IS] = sort(roundn10(d,round(log10(epsilon))+3));
		d = d(IS).';
		Q(:,j:k) = Q(:,j:k)*V(:,IS);
		IE(j:k) = [abs(d(1:end-1)-d(2:end)) < epsilon,0];
	end
end

%%
% Replace the 'roundn' function from the mapping toolbox
% to avoid package dependency
function x = roundn10(x,n)

factors = 10^(fix(-n));
x = round(x*factors)/factors;
