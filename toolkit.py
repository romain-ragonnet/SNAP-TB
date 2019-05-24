from numpy import exp
from pyDOE import lhs

def make_sigmoidal_curve(y_low=0, y_high=1., x_start=0, x_inflect=0.5, multiplier=1.):
    """
    Base function to make sigmoidal curves.

    Args:
        y_low: lowest y value
        y_high: highest y value
        x_inflect: inflection point of graph along the x-axis
        multiplier: if 1, slope at x_inflect goes to (0, y_low), larger
                    values makes it steeper
    Returns:
        function that increases sigmoidally from 0 y_low to y_high
        the halfway point is at x_inflect on the x-axis and the slope
        at x_inflect goes to (0, y_low) if the multiplier is 1.
    """

    amplitude = y_high - y_low
    if amplitude == 0:
        def curve(x):
            return y_low
        return curve

    x_delta = x_inflect - x_start
    slope_at_inflection = multiplier * 0.5 * amplitude / x_delta
    b = 4. * slope_at_inflection / amplitude

    def curve(x):
        arg = b * ( x_inflect - x )
        # check for large values that will blow out exp
        if arg > 10.0:
            return y_low
        return amplitude / ( 1. + exp( arg ) ) + y_low

    return curve

def make_scale_up_function(x,y):


    def fn(t):
        if t < 100.:
            y = 0.
        elif t > 200.:
            y = 1.
        else:
            y = (t-100.)/100.

        return y

    return fn

def lhs_sampler(n_params, n_samples):
    out = lhs(n=n_params, samples=n_samples, criterion='c')
    return out
