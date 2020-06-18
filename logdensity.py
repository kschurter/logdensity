"""Estimate a local polynomial approximation to a density

This module contains just two functions that implement the estimator
developed in Pinkse, J. and Schurter, K. (2020) "Estimates of
derivatives of (log) densities and related objects." Available from
<http://arxiv.org/abs/2006.01328>.

Version: 0.1.0 (Jun 18 2020)
Depends:
    numpy
    scipy
License:  GNU General Public License version 3
Contact:  Karl Schurter <kschurter@psu.edu>
Examples:
    import logdensity as ld
    import random

    gaussdat = [random.gauss(0,1) for i in range(1000)]  # iid data
    x = [0,0.1,0.2]  # evaluation points
    h = 0.5  # bandwidth
    ld.logdensity(gaussdat,x,h)
    ld.logdensity(gaussdat,x,h,S=2)

    expodat = [random.expovariate(1) for i in range(1000)]
    ld.logdensity(expodat,x,h,minx=0)
"""

import warnings
import traceback

import numpy as np
import scipy.linalg as linalg
import scipy.integrate as integrate
from scipy.special import factorial

warnings.filterwarnings("error",".",linalg.LinAlgWarning)
_msg = ('\nWarnings filter modified to raise LinAlgWarning as error.\n'
        + 'This is recommended for use with logdensity() and logdensity_fit()\n'
        + 'because scipy.linalg._solve_check() sometimes fails to detect\n'
        + 'errors caused by bad user input. Revert to default at your own '
        + 'risk.')
warnings.warn(_msg,Warning)

def logdensity_fit(data,x,h,g,dg,minx=float('-inf'),maxx=float('inf'),logf=True,
                   m=None,quadrature_args={},exact=True,*args,**kwargs):
    """The workhorse for logdensity().

    This function has fewer defaults than logdensity(), does not do any
    argument verification, or accommodate vector-valued x or h. This
    function is intended to be used if you have read the documentation for
    logdensity(), are confident your g and dg functions meet the necessary
    criteria, and you now care about performance.

    Args:
        logdensity_fit() accepts the same arguments that logdensity()
        accepts, except for S, which is inferred from the shape of the
        output of g. See help(logdensity) for further explanation of the
        following arguments, but defer to the types documented here:

        data:  a 1-D numpy array. nan is not allowed.
        x:  a finite float.
        h:  a positive, finite float.
        g:  a function.
        dg:  a function.
        minx:  a float, possibly -inf.
        maxx:  a float, possibly inf.
        logf:  logical.
        m:  a function of a single vector argument that returns a list-like
            value the same size as its argument.
        quadrature_args:  a dictionary.

        logdensity_fit() also accepts

        exact: logical. Specifies whether an exact solution to the integral
            in the estimate of the log-density should be used (if possible).
            Currently only works with local linear estimation (S=1) and the
            Epanechnikov kernel ( m(u) = 0.75*(1-u**2)*(abs(u)<=1) ).

    Returns:
        A 1-D numpy array containing the estimated log-density and its
        derivative(s). The j-th index is the estimated j-th derivative of
        the natural logarithm of the density of the data. The 0-th index is
        numpy.nan if not logf.

    Raises:
        ValueError
            If g or dg returns a non-finite value.
        LinAlgError
            If a singular matrix is passed to scipy.linalg.solve(). This
            might mean there are not enough observations within one
            bandwidth to estimate as many derivatives as you are requesting. 
        LinAlgWarning
            If an ill-conditioned matrix is by scipy.linalg.solve(). On
            import, this module changes the warnings filter to raise this
            warning as an error.
        RuntimeWarning
            Any exceptions thrown by scipy.integrate.quadrature() will be
            caught and reissued as a warning. This is fatal to the estimates
            of the log-density, only. The estimates of the derivatives are
            unaffected, so the function will continue.           
    """

    u = (data-x)/h
    zl = min((x-minx)/h, 1)
    zr = min((maxx-x)/h, 1)
    u = u[( u>= -zl) & (u <= zr)]
    gu = g(u,zl,zr,*args,**kwargs)
    rangeS = np.array(range(gu.shape[0]))
    R = -1 * np.dot(gu, u[:,None]**rangeS)
    el = np.sum(dg(u,zl,zr,*args,**kwargs), axis = 1)
    try:
        beta = linalg.solve(R,el) * factorial(rangeS) / (h**(rangeS+1))
    except (linalg.LinAlgError, linalg.LinAlgWarning) as exc:
        moreinfo = ('Attempted to estimate %d derivative(s) at %g '%(el.size,x)
                    + 'using %d observations. If that sounds '%(u.size)
                    + 'impossible, consider increasing h from %g.'%h,)
        exc.args += moreinfo
        raise
    except ValueError as vexc:
        moreinfo = ()
        if np.any(~np.isfinite(R)):
            moreinfo += ('Message from logdensity_fit(): Encountered '
                         + 'non-finite output from g(u,%g,%g).'%(zl,zr),)
        if np.any(~np.isfinite(el)):
            moreinfo += ('Message from logdensity_fit(): Encountered '
                         + 'non-finite output from dg(u,%g,%g).'%(zl,zr),)
        vexc.args +=moreinfo
        raise
    if logf:
        if exact and (beta.size==1):
            num = sum((1-u**2)) / len(data)
            def denfun(x,b):
                return np.exp(b*x) * (b**2 * (1-x**2) + 2*b*x - 2) / b**3
            den = denfun(zr,beta*h) - denfun(-zl,beta*h)
        else:
            num = sum(m(u)) / len(data)
            def denint(x):
                return m(x)*np.exp(np.dot((x[:,None]*h)**(rangeS+1),
                                          beta/factorial(rangeS+1)))
            try:
                den = integrate.quadrature(denint,-zl,zr,**quadrature_args)[0]
            except Exception as exc:
                msg = ('Numerical integration failed when trying to estimate '
                       + 'the log-density at %g. '%x
                       + 'Derivative estimates are unaffected.\n'
                       + '%s'%traceback.format_exc(1))
                warnings.warn(msg,RuntimeWarning)
                den = np.nan
        ld = np.log(num) - np.log(den) - np.log(h)
    else:
        ld = np.nan
    return np.append(ld,beta)


def logdensity(data,x,h,g=None,dg=None,S=1,minx=float('-inf'),maxx=float('inf'),
               logf=True,m=None,quadrature_args={},*args,**kwargs):
    """ Local polynomial estimation of the log density

    Estimates the (natural) logarithm of a density and its derivatives from
    iid data using a local polynomial approximation. The function is based
    on Pinkse, J. and Schurter, K. (2020) "Estimates of derivatives of (log)
    densities and related object." The paper is available on arXiv:
    <http://arxiv.org/abs/2006.01328>.

    Args:
        data: a list-like object containing iid observations. Will be made
            into a 1-D numpy array using numpy.ravel().
        x:  a list-like object containing the points at which to estimate
            the log-density and its derviative(s). x will also be ravel'ed.
        h:  a list-like object containing the bandwidth(s) to be used in the
            estimation at the corresponding point x. h may have length 1,
            in which case the same "fixed" bandwidth is used for all x. If
            the length of x is a multiple of the length of h, h will be
            repeated to match. h will also be ravel'ed before using.
        g:  a function (defined in the paper) that will be used to estimate
            the derivatives of the log-density. The function must accept at
            least three positional arguments
                1. u: a vector of points at which to evaluate the function,
                2. zl: a scalar specifying the left end of g's support,
                3. zr: a scalar specifying the right end of g's support.
            The number of columns of the output of g must match the size of
            its first argument. If the output has S rows, a S-degree local
            polynomial approximation to the log-density will be estimated. g
            function must equal numpy.zeros(S) at the left and right ends of
            its support [-zl,zr] contained in [-1,1]. The j-th entry of the
            default g(u,zl,zr) equals (u+zl)**j * (u-zr).
        dg: a function the returns the derivative of g with respect to its
            first argument.
        S:  a positive integer specifying the degree of the local polynomial
            approximation. This argument only matters if users do not supply
            their own g and dg function, otherwise it is ignored.
        minx:  the left end of the support of the data. Defaults to '-inf'.
        maxx:  the right end of the support of the data. Defaults to 'inf'.
        logf:  logical. compute the logarithm of the density or just its
            derivative(s).
        m:  a kernel function to use to estimate the log-density. m should
            accept one argument and return a value of the same size as that
            argument. Defaults to m(u) = (1-u**^2)*(abs(u)<=1).
        quad_args:  a dictionary of keyword arguments passed to
            scipy.integrate.quadrature()
        *args:  further positional arguments that passed to g and dg.
        **kwargs:  further keyword arguments that passed to g and dg.

    Returns:
        A numpy array containing the estimated log-density (if logf) and its
        derivatives. The [,j] column of the array contains the estimates of
        the j-th derivative of the log-density evaluated at the
        corresponding element of numpy.ravel(x).

    Raises:
        RuntimeWarning
        UserWarning
        LinAlgWarning
        LinAlgError
        ValueError           
    """

    # Validate arguments and coerce data, x, and h into 1-D numpy arrays.
    data = np.ravel(data)
    x = np.ravel(x)
    h = np.ravel(h)
    if any(~np.isfinite(h)) | any(h<=0):
        raise ValueError('Bandwidth h must be positive and finite.')
    if any(np.isnan(data)):
        nans = np.isnan(data)
        data = data[~nans]
        msg = '%d missing observation(s) dropped.'%sum(nans)
        warnings.warn(msg)
    if any(data > maxx) or any(data < minx):
        raise ValueError('All data must be in the support [%g,%g].'%(minx,maxx))
    if any(x < minx) or any(x > maxx) or not all(np.isfinite(x)):
        msg = ('Ignoring infinite values x and values of x that lie outside '
              + 'the support [%g,%g]'%(minx,maxx))
        warnings.warn(msg)
        x = x[(x>=minx) & (x<=maxx) & np.isfinite(x)]
        if x.size==0:
            raise ValueError('No finite values of x within support' +
                            '[%g,%g]'%(minx,maxx))

    # Define default g, dg, and m. If user supplied g and dg, check
    # restrictions on g.
    if (g is None) or (dg is None):
        if (g is not None) or (dg is not None):
            raise ValueError('User must specify both g and dg or neither.')
        if ((S % 1) != 0) or (S<1):
            raise ValueError('S must be a positive integer.')
        # g and dg will not be evaluated at points u outside [-zl,zr],
        # So we're not worrying about making sure they evaluate to zero
        # on (-inf,-zl) and (zr,inf).
        def g(u,zl,zr,S):
            s = np.expand_dims(range(1,S+1), axis = -1)
            return ((u+zl)**s) * (u-zr)
        def dg(u,zl,zr,S):
            s = np.expand_dims(range(S), axis = -1)
            return ((u+zl)**s) * ((u+zl) + ((u-zr)*(s+1)))
        args = [S]
    else:
        zl = np.minimum((x-minx)/h, 1)
        zr = np.minimum((maxx-x)/h, 1)
        gatboundaries = [g([-zl,zr],zl,zr) for zl,zr in zip(zl,zr)]
        if np.any(np.array(gatboundaries)!=0):
            raise ValueError('g must be zero at the boundary of its support.')
    exact = (m is None)
    if exact:
        def m(u):
            return (1-u**2)*(abs(u)<=1)

    # Apply logdensity_fit to pairs of (x, h)
    if h.size < x.size:
        h = np.repeat(h,x.size//h.size)
    arr = [logdensity_fit(data,x_el,h_el,g,dg,minx,maxx,logf,m,quadrature_args,
                          exact,*args,**kwargs)
           for x_el, h_el in zip(x,h)]
    return np.array(arr)
