from numpy import array

# Simple parameters storage
class Parameters:
    def __init__(s, q0: array):
        # The initial configration
        s.q0 = q0
        # Gravity vector
        s.g = array([0, 0, 9.81])
        # Friction coefficient
        s.mu = 0.7
        # Position error gain
        s.Kp = 200
        # Max swing height
        s.swing_height = 0.15
        # Frequency to interpolate the solution to
        s.desired_frequency = 0.01
        # Degree of interpolating polynomial
        s.degree = 1
        # Problem options for casadi
        s.p_opts = {"expand": True}
        # Solver options for casadi
        s.s_opts = {"max_iter": 10000, 
                    # "linear_solver": "ma97"
                    }