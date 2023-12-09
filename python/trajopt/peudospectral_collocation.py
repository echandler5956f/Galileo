import numpy as np
import casadi as ca


# Class for setting up pseudospectral collocation
class PseudoSpectralCollocation:
    def __init__(s, degree: int):
        # Degree of interpolating polynomial
        s.degree = degree

        # Get collocation points
        tau_root = np.append(0, ca.collocation_points(s.degree, "radau"))

        # Coefficients of the collocation equation
        s.C = np.zeros((s.degree + 1, s.degree + 1))

        # Coefficients of the continuity equation
        s.D = np.zeros(s.degree + 1)

        # Coefficients of the quadrature function
        s.B = np.zeros(s.degree + 1)

        # Construct polynomial basis
        for j in range(s.degree + 1):
            # Construct Lagrange polynomials to get the polynomial basis at the collocation point
            p = np.poly1d([1])
            for r in range(s.degree + 1):
                if r != j:
                    p *= np.poly1d([1, -tau_root[r]]) / (tau_root[j] - tau_root[r])

            # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            s.D[j] = p(1.0)

            # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the
            # continuity equation
            p_der = np.polyder(p)
            for r in range(s.degree + 1):
                s.C[j, r] = p_der(tau_root[r])

            # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
            pint = np.polyint(p)
            s.B[j] = pint(1.0)
