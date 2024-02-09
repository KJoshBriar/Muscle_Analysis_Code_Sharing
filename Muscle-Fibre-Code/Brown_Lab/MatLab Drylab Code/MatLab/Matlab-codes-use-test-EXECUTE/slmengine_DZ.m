function slm = slmengine_DZ(x,y,varargin)
% This function was is the same as the slmengine, however I went through it
% and I removed all the extra settings that I am not using. The reason I
% did this was so that I could run it faster and more optimized for what
% I'm doing. It gives the same answer as slmengine with 10x the speed. The
% only adjustable parameters are the plot ('on' or 'off'(default) and the
% number of equally spaced knots (defaut is 6). 

% Edited by Derek Zwambag on April 15, 2018

%% Parse input parameters
prescription.knots = 6;
prescription.plot = 'off';

if ~isempty(varargin)
    nparam = length(varargin)/2;
    if nparam ~= floor(nparam)
        error('Parameters must come in pairs')
    end
    
    for k = 1:nparam
        if strcmp(varargin{k*2-1},'knots')
            prescription.knots = varargin{k*2};
        elseif strcmp(varargin{k*2-1},'plot')
            prescription.plot = varargin{k*2};
        else
            error('%s is not a valid input',varargin{k*2-1})
        end
    end
end

%% check the data for size, turning it into column vectors
x = x(:);
y = y(:);
n = length(x);
if n~=length(y)
  error('SLMENGINE:inconsistentdata','x and y must be the same size')
end

%% were there any NaN or inf elements in the data?
k = isnan(x) | isnan(y) | isinf(x) | isinf(y);
if any(k)
  % drop them from the analysis
  x(k) = [];
  y(k) = [];
end

% we need to scale y to minimize any numerical issues.
% note that the scaling does not change the signs of any derivatives,
% so monotonicity and curvature constraints are not accidentally
% inverted.
% scaleproblem modifies y into yhat, so that yhat now lies in the
% interval [1/phi,phi], where phi is the golden ratio. This
% transformation has the property that all dependent values vary
% over a range of 1.
[yhat,prescription,YScale,YShift] = scaleproblem(y,prescription);

% The knots are fixed. Estimate the model
slm = slmengine_cubic(x,yhat,prescription);

% back out any shift & scale factors from the model
slm.coef(:,1) = (slm.coef(:,1) - YShift)/YScale;
if size(slm.coef,2) > 1
    slm.coef(:,2) = slm.coef(:,2)/YScale;
end

slm.x = x;
slm.y = y;

% do we need to plot the curve? I.e., what was the
% prescription.Plot setting?
if strcmp(prescription.plot,'on')
    % plot the curve
    plotslm(slm)
end

end % End of main function slmengine


% ========================================================
% ========= set scaling as necessary =========
% ========================================================
function [yhat,prescription,YScale,YShift] = scaleproblem(y,prescription)

% chooses an appropriate scaling for the problem
% there is no need to scale x, since each knot interval is already
% implicitly normalized into [0,1]. x is passed in only for
% a few parameters, such as the integral constraint, where x would
% be needed.

% scale y so that the minimum value is 1/phi, and the maximum value phi.
% where phi is the golden ratio, so phi = (sqrt(5) + 1)/2 = 1.6180...
% Note that phi - 1/phi = 1, so the new range of yhat is 1. (Note that
% this interval was carefully chosen to bring y as close to 1 as
% possible, with an interval length of 1.)
%
% The transformation is:
%   yhat = y*yscale + yshift
phi_inverse = (sqrt(5) - 1)/2;

% shift and scale are determined from the min and max of y.
ymin = min(y);
ymax = max(y);

YScale = 1./(ymax - ymin);
if isinf(YScale)
    % in case data was passed in that is constant, then
    % the range of y is zero. No scaling need be done then.
    YScale = 1;
end

% recover the shift once the scale factor is known.
YShift = phi_inverse - YScale*ymin;

% scale y to refect the shift and scale
yhat = y*YScale + YShift;
  
prescription.YScale = YScale;
prescription.YShift = YShift;
end

% ========================================================
% =============== piecewise cubic model ==================
% ========================================================
function slm = slmengine_cubic(x,y,prescription)
% fits a piecewise cubic shape prescriptive model

% slmengine has already made 'x' a column vector, and ensured compatibility
% between the length of 'x' and 'y'
nx = length(x);

% knots vector, dx
if length(prescription.knots) == 1
    nk = prescription.knots;
    knots = linspace(min(x),max(x),nk)';
    knots(end) = max(x); % just to make sure
else
    knots = sort(prescription.knots(:));
    nk = length(knots);
    if (knots(1) > min(x)) || (knots(end) < max(x))
        error('SLMENGINE:inadequateknots',['Knots do not contain the data. Data range: ',num2str([min(x),max(x)])])
    end
end
dx = diff(knots);
if any(dx==0)
    error('SLMENGINE:indistinctknots','Knots must be distinct')
end

% number of coefficients to estimate.
% a piecewise cubic Hermite has two coefficients at each knot.
nc = 2*nk;

% create empty equality constraints arrays and rhs vectors
Meq = zeros(0,nc);
rhseq = [];

% -------------------------------------
% build design matrix - 
% first, bin the data - histc will do it best
[junk,xbin] = histc(x,knots); %#ok
% any point which falls at the top end, is said to
% be in the last bin.
xbin(xbin==nk)=nk-1;

% build design matrix
t = (x - knots(xbin))./dx(xbin);
t2 = t.^2;
t3 = t.^3;
s2 = (1-t).^2;
s3 = (1-t).^3;

vals = [3*s2-2*s3 ; 3*t2-2*t3 ; ...
  -(s3-s2).*dx(xbin) ; (t3-t2).*dx(xbin)];

% the coefficients will be stored in two blocks,
% first nk function values, then nk derivatives.
Mdes = accumarray([repmat((1:nx)',4,1), ...
  [xbin;xbin+1;nk+xbin;nk+xbin+1]],vals,[nx,nc]);
rhs = y;

% -------------------------------------
% Regularizer
% We are integrating the piecewise linear f''(x), as
% a quadratic form in terms of the (unknown) second
% derivatives at the knots.
Mreg=zeros(nk,nk);
Mreg(1,1:2)=[dx(1)/3 , dx(1)/6];
Mreg(nk,nk+[-1 0])=[dx(end)/6 , dx(end)/3];
for i=2:(nk-1)
  Mreg(i,i+[-1 0 1])=[dx(i-1)/6 , (dx(i-1)+dx(i))/3 , dx(i)/6];
end
% do a matrix square root. cholesky is simplest & ok since regmat is
% positive definite. this way we can write the quadratic form as:
%    s'*r'*r*s,
% where s is the vector of second derivatives at the knots.
Mreg=chol(Mreg);

% next, write the second derivatives as a function of the
% function values and first derivatives.   s = [sf,sd]*[f;d]
sf = zeros(nk,nk);
sd = zeros(nk,nk);
for i = 1:(nk-1)
  sf(i,i+[0 1]) = [-1 1].*(6/(dx(i)^2));
  sd(i,i+[0 1]) = [-4 -2]./dx(i);
end
sf(nk,nk+[-1 0]) = [1 -1].*(6/(dx(end)^2));
sd(nk,nk+[-1 0]) = [2 4]./dx(end);
Mreg = Mreg*[sf,sd];

% scale the regularizer before we apply the
% regularization parameter.
Mreg = Mreg/norm(Mreg,1);
rhsreg = zeros(nk,1);

% -------------------------------------
% C2 continuity across knots
MC2 = zeros(nk-2,nc);
for i = 1:(nk-2)
    MC2(i,[i,i+1]) = [6 -6]./(dx(i).^2);
    MC2(i,[i+1,i+2]) = MC2(i,[i+1,i+2]) + [6 -6]./(dx(i+1).^2);

    MC2(i,nk+[i,i+1]) = [2 4]./dx(i);
    MC2(i,nk+[i+1,i+2]) = MC2(i,nk+[i+1,i+2]) + [4 2]./dx(i+1);
end

Meq = [Meq;MC2];
rhseq = [rhseq;zeros(nk-2,1)];

% -------------------------------------
% scale equalities for unit absolute row sum
rs = diag(1./sum(abs(Meq),2));
Meq = rs*Meq;
rhseq = rs*rhseq;

% -------------------------------------
% now worry about the regularization. 
% 1. We have a given regularization parameter

% Regularization parameter
RP = 0.0001; % default value

% combine design matrix and regularizer
Mfit = [Mdes;RP*Mreg];
rhsfit = [rhs;RP*rhsreg];

% with no inequality constraints, lse is faster than
% is lsqlin. This also allows the use of slm when
% the optimization toolbox is not present if there
% are no inequality constraints.

A = Mfit; b = rhsfit; C = full(Meq); d = rhseq;

[n,p] = size(A);
[r,nrhs] = size(b);
[m,ccols] = size(C);
if n~=r
  error 'A and b are incompatible in size (wrong number of rows)'
elseif ~isempty(C) && (p~=ccols)
  error 'A and C must have the same number of columns'
elseif ~isempty(C) && issparse(C)
  error 'C may not be a sparse matrix'
elseif ~isempty(C) && (m~=size(d,1))
  error 'C and d are incompatible in size (wrong number of rows)'
elseif ~isempty(C) && (size(d,2)~=1)
  error 'd must have only one column'
elseif isempty(C) && ~isempty(d)
  error 'C and d are inconsistent with each other (one was empty)'
elseif ~isempty(C) && isempty(d)
  error 'C and d are inconsistent with each other (one was empty)'
end


% tolerance used on the solve
Ctol = 1.e-13;

% allow a rank deficient equality constraint matrix
% column pivoted qr to eliminate variables
[Q,R,E]=qr(C,0);

% get the numerical rank of R (and therefore C)
if m == 1
%      rdiag = R(1,1);
  rdiag = abs(R(1,1));
else
  rdiag = abs(diag(R));
end
crank = sum((rdiag(1)*Ctol) <= rdiag);
if crank >= p
  error 'Overly constrained problem.'
end

% check for consistency in the constraints in
% the event of rank deficiency in the constraint
% system
if crank < m
  k = Q(:,(crank+1):end)'*d;
  if any(k > (Ctol*norm(d)));
    error 'The constraint system is deficient and numerically inconsistent'
  end
end

% only need the first crank columns of Q
qpd = Q(:,1:crank)'*d;

% which columns of A (variables) will we eliminate?
j_subs = E(1:crank);
% those that remain will be estimated
j_est = E((crank+1):p);

r1 = R(1:crank,1:crank);
r2 = R(1:crank,(crank+1):p);

A1 = A(:,j_subs);
qpd = qpd(1:crank,:);

% note that \ is still ok here, even if pinv
% is used for the main regression.
bmod = b-A1*(r1\repmat(qpd,1,nrhs));
Amod = A(:,j_est)-A1*(r1\r2);

% now solve the reduced problem without weights
x2 = Amod\bmod;

% recover eliminated unknowns
x1 = r1\(repmat(qpd,1,nrhs)-r2*x2);

% stuff all estimated parameters into final vector
coef = zeros(p,nrhs);
coef(j_est,:) = x2;
coef(j_subs,:) = x1;   

% -------------------------------------
% unpack coefficients into the result structure
slm.form = 'slm';
slm.degree = 3;
slm.knots = knots;
slm.coef = reshape(coef,nk,2);
end