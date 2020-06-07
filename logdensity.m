function logden = logdensity(data,x,h,varargin)
% LOGDENSITY  Estimate the log-density and its derivatives
%
%   LOGDENSITY(data,x,h) estimates the logarithm of the density of the iid
%   data at x using bandwidth h. The parameters data and x are numeric
%   arrays that will be coerced into vectors. The bandwidth must be a
%   positive, finite scalar. The estimates are returned in a matrix with
%   length(x) columns and S+1 rows (default S = 1).
%
%   LOGDENSITY(data,x,h,PARAM1,VAL1,PARAM2,VAL2,...) estimates the
%   log-density and derivatives with optional parameters. These options are
%
%   'g',  a function
%   'dg', a function that returns the derivatve of 'g'.
%         'g' and 'dg' are used to estimate the derivatives of the 
%         log-density. The functions must accept four arguments:
%             1. a vector of points at which to evaluate the function,
%             2. the left end of the support of the function 'g'
%             3. the right end of the support of the funciton 'g'
%             4. the order of the polynomial approximation 'S'.
%         'g' must equal zeros(S,1) at the boundaries of its support
%         [-zl,zr] which is contained in [-1,1]. 'dg' is the derivative of
%         'g' with respect to its first argument. The default g(u,zl,zr) is
%         a function that returns an S-length vector whose j-th element is
%         (u + zl).^j .* (zr-u), and 'dg' defaults to its derivative with 
%         respect to u.
%
%   'S', a positive integer indicating the order of the local polynomial
%        approximation to the logarithm of the density
%
%   'minx', the left end of the support of the density
%
%   'maxx', the right end of the support of the density
%   
%   'mz', a scalar-valued function of a single variable
%         'mz' is used to estimate the log-density. The default value is
%         the epanechnikov kernel mz(u) = 0.75*(1-u.^2) .* (abs(u) <= 1).
%
%   'logf', logical indicating whether to estimate the density itself or
%           just return the derivative estimates. If false, the first row
%           of the result will be NaN.
%
%   Example:
%      dat = chi2rnd(2,100)
%      logdensity(dat,[0, 0.1, 0.2],0.5,'minx',0)
%
%   References:
%       This function is based on Pinkse, J. and Schurter, K. (2020)
%       "Estimates of derivatives of (log) densities and related objects."
%       available at <https://arxiv.org/abs/2006.01328>.
  
  %% Parse arguments and set defaults for optional parameters
  
  p=inputParser;
  addRequired(p,'data');               % data points
  addRequired(p,'x');                  % evaluation points
  addRequired(p,'h');                  % bandwidth
  addParameter(p,'S',1);               % order of polynomial approximation
  addParameter(p,'g',@(u,zl,zr,S) bsxfun(@times,bsxfun(@power,(u+zl), (1:S)),(zr-u)));
  addParameter(p,'dg',@(u,zl,zr,S) bsxfun(@power,(u+zl),(0:(S-1))) .* bsxfun(@minus, bsxfun(@times,(zr-u),1:S),(u + zl)));
  addParameter(p,'minx',-Inf);         % bottom of support
  addParameter(p,'maxx',Inf);          % top of support
  addParameter(p,'mz',@(x) (abs(x)<=1).*0.75.*(1-x.*x));   % kernel for estimating log-density
  addParameter(p,'logf',true);         % estimate the log-density, too?
  
  parse(p,data,x,h,varargin{:});
  dat=p.Results.data(:);
  x=p.Results.x(:);
  h=p.Results.h;
  g=p.Results.g;
  dg=p.Results.dg;
  S=p.Results.S;
  minx=p.Results.minx;
  maxx=p.Results.maxx;
  mz=p.Results.mz;
  logf=p.Results.logf;
  
  %% Do some initial parameter validation
  if (S<1 || (floor(S)~=S))
      error('Polynomial order S must be a positive integer.')
  elseif (0>h || ~isfinite(h))
      error('Bandwidth must be a positive finite number.')
  elseif any(isnan(dat))
      drop = isnan(dat);
      dat = dat(~drop);
      warning('Dropped %d missing observation(s)',sum(drop))
  elseif any(minx>dat | dat>maxx)
      error('Not all data points are in the support [minx,maxx].')
  elseif any(minx>x | x >maxx | ~isfinite(x)) 
      warning('Some evaluation point(s) are infinite or not in the support [minx,maxx]. Ignoring offending elements of x.')
      x = x(x>=minx & x<= maxx & isfinite(x));
      if isempty(x)
          error('No finite evaluation points are in the support.')
      end
  elseif xor(any(strcmp('g',p.UsingDefaults)),any(strcmp('dg',p.UsingDefaults)))
      error('User must supply both g and dg or neither.')
  end
  
  %% Compute estimates
  nx=length(x);
  ndat=length(dat);
  zl=min((x-minx)/h,1);
  zr=min((maxx-x)/h,1);
  logden=NaN(S+1,nx);
  v = bsxfun(@minus, dat, x') ./ h;
  for ix=1:nx
    u = v(v(:,ix)>-zl(ix) & v(:,ix)<zr(ix),ix);
    if (S > length(u))
        error('Not enough observations within a bandwidth of x to estimate S derivatives.')
    elseif any(g(-zl(ix),zl(ix),zr(ix),S)~=0)
       error('Function g used is not zero at lower end of its support.')
    elseif any(g(zr(ix),zl(ix),zr(ix),S)~=0)
       error('Function g used is not zero at upper end of its support.')
    end       
    l=sum(dg(u,zl(ix),zr(ix),S), 1);
    if any(~isfinite(l))
        error('Function dg returned a non-finite value')
    end
    R=(g(u,zl(ix),zr(ix),S)') * bsxfun(@power,u,0:(S-1));
    if any(~isfinite(R),'all')
        error('Function g returned a non-finite value')
    end
    logden(2:end,ix)=-(R\(l')).*(factorial(0:(S-1))./(h.^(1:S)))';
    if logf 
       num=sum(mz(u))/(ndat*h); % numerator of estimate of density
       if S==1 && any(strcmp('mz',p.UsingDefaults)) % exact solution for denominator
           denfunc = @(x,c) (0.75*exp(c*x)*(c*c*(1-x*x)+2*c*x-2)/c^3);  
           denom=denfunc(zr(ix),logden(2,ix)*h)-denfunc(-zl(ix),logden(2,ix)*h);
       else % numerical solution for denominator
           dint = @(u) (mz(u).* exp((logden(2:end,ix)'./factorial(1:S)) * bsxfun(@power, (u.*h), (1:S)')));
           try
              denom=quadgk(dint,-zl(ix),zr(ix));
           catch err
              denom = NaN;
              warning('Numerical integration failure in estimate of the log-density at %0.4g. Derivative estimates are unaffected.',x(ix))
              warning('Information from quadgk():\n identifier =  %s\n  message = %s',err.identifier,err.message) 
           end
       end
       logden(1,ix) = log(num) - log(denom);
    end   
  end
end