import numpy as np
import casadi as ca

class ModelSymbols:
    def __init__(s, nq: int, nv: int):
        # Size of joint position vector
        s.nq = nq
        # Size of joint tangent vector
        s.nv = nv
        # 3 linear, 3 angular
        s.nm = 3
        # Size of centroidal momentum vector
        s.nh = 2 * s.nm
        # Size of centroidal momentum time derivative vector
        s.ndh = 2 * s.nh

        # Size of the state vector
        s.nx = 2 * s.nh + s.nq + s.nv
        # Size of the state delta vector
        s.ndx = 2 * s.nh + 2 * s.nv
        # State vector where x = [h; hdot; q; v]
        s.cx = ca.SX.sym("x", s.nx, 1)
        # State delta vector
        s.cdx = ca.SX.sym("dx", s.ndx, 1)

        # Contact force and contact wrench
        s.cu = ca.SX.sym("u", 6, 1)
        # Joint tangents are one of the decision variables
        s.cvju = ca.SX.sym("vju", s.nv - 6, 1)

        # Momenta: nh x 1
        s.ch = s.cx[: s.nh]
        # Momenta delta: nh x 1
        s.ch_d = s.cdx[: s.nh]

        # Momenta time derivative: nh x 1
        s.cdh = s.cx[s.nh : s.ndh]
        # Momentum time derivative delta: nh x 1
        s.cdh_d = s.cdx[s.nh : s.ndh]

        # q: nq x 1
        s.cq = s.cx[s.ndh : s.nq + s.ndh]
        # q delta: nv x 1
        s.cq_d = s.cdx[s.ndh : s.nv + s.ndh]

        # qj: (nq - 7) x 1
        s.cqj = s.cq[7:]

        # v: nv x 1
        s.cv = s.cx[s.ndh + s.nq :]
        # v delta: nv x 1
        s.cv_d = s.cdx[s.ndh + s.nv :]

        # v_j: (nv - 6) x 1
        s.cvj = s.cv[6:]

        # f: 3 x 1
        s.cf = s.cu[:3]
        # tau: 3 x 1
        s.ctau = s.cu[3:]

        # Polynomial variables
        s.cxp = ca.SX.sym("xp", s.nx, 1)
        # Polynomial momenta: nh x 1
        s.chp = s.cxp[: s.nh]
        # Polynomial momentum time derivative: nh x 1
        s.cdhp = s.cxp[s.nh : s.ndh]
        # Polynomial q: nq x 1
        s.cqp = s.cxp[s.ndh : s.ndh + s.nq]
        # Polynomial v: nv x 1
        s.cvp = s.cxp[s.ndh + s.nq :]

        # 3D vector helper
        s.c3vec = ca.SX.sym("3vec", 3, 1)
        # Quaternion vector helper
        s.cqvec = ca.SX.sym("qvec", 4, 1)

        # Time helper variable
        s.cdt = ca.SX.sym("dt", 1, 1)