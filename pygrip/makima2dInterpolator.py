'''
This is the Makima2DInterpolator class.
It allows to define an interpolator based on the Modified Akima method.
A linear extrapolation method is used across the bounds.
'''
import numpy as np
from scipy.interpolate import PchipInterpolator, Akima1DInterpolator

class Makima2DInterpolator:

    # Constructor
    def __init__(self, x, y, z):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.z = np.asarray(z)

        if self.z.shape != (len(y), len(x)):
            raise ValueError(f"Z shape must be ({len(y), len(x)}), instead got {self.z.shape}")

    
    def _akima_with_extrapolation(self, x, y_vals):
        # Return an interpolator with linear extrapolaion beyond bounds
        akima = Akima1DInterpolator(x, y_vals, extrapolate=True)
        def extrap(xq):
           xq = np.asarray(xq)
           result = akima(xq)
           nfit = 16
        #    coeffs = np.polyfit(x, y_vals, deg=3)
           # Left extrapolation
           left = xq < x[0]
           if np.any(left):
               coeffs = np.polyfit(x[:nfit], y_vals[:nfit], deg=2)
               result[left] = np.polyval(coeffs, xq[left])
                # slope = (y_vals[1] - y_vals[0]) / (x[1] - x[0])
                # result[left] = y_vals[0] + slope * (xq[left] - x[0])
            
            # Right extrapolation
           right = xq > x[-1]
           if np.any(right):
               coeffs = np.polyfit(x[-nfit:], y_vals[-nfit:], deg=2)
               result[right] = np.polyval(coeffs, x[right])
            #    slope = (y_vals[-1] - y_vals[-2]) / (x[-1] - x[-2])
            #    result[right] = y_vals[-1] + slope * (xq[right] - x[-1])
            
           return result
        return extrap
    
    # Call
    def __call__(self, xq, yq):

        # Interpolate along x (axis=1) for each value of y
        intermediate = []
        for row in self.z:
            interp = self._akima_with_extrapolation(self.x, row)
            intermediate.append(interp(xq))
        intermediate = np.array(intermediate)

        # Interpolate along y (axis=0)
        interp_final = self._akima_with_extrapolation(self.y, intermediate)
        return float(interp_final(yq))