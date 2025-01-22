import numpy as np
import sympy as sp
from Node import No
from math import sqrt, sin, cos


class Elemento:
    """
    Class that implements a finite element for a frame structure.
    Assumes linearly distributed loads (q1 to q2 and p1 to p2).
    """

    def __init__(self, A: No, B: No, EA, EI, Id, q1=0, q2=0, p1=0, p2=0):
        # assume linear loads
        self.A = A
        self.B = B
        self.Id = Id
        # Angle of element with x-axis and length
        self.x = B.x - A.x
        self.y = B.y - A.y
        self.L = sqrt(self.x ** 2 + self.y ** 2)
        self.alfa = np.arctan2(self.y, self.x)
        self.EA = EA
        self.EI = EI
        self.p1 = p1
        self.p2 = p2
        self.q1 = q1
        self.q2 = q2
        # rotation matrix
        c = cos(self.alfa)
        s = sin(self.alfa)
        self.R = np.array([
            [c, -s, 0, 0, 0, 0],
            [s, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, -s, 0],
            [0, 0, 0, s, c, 0],
            [0, 0, 0, 0, 0, 1]
        ])
        # Normal force, shear force and bending moment.
        self.N = None
        self.V = None
        self.M = None

        # Stiffness matrix and forces
        self.f = np.zeros((6, 1))
        self.K = np.zeros((6, 6))
        self.calcularMatrizdeRigidez()
        self.calcularMatrizdeForcas()

    def calcularMatrizdeRigidez(self):
        """
        Computes Stiffness matrix
        """
        K_trelica = np.zeros((6, 6))
        for i in [0, 3]:
            for j in [0, 3]:
                K_trelica[i, j] = self.EA / self.L * (-1) ** (i + j)
        K_viga = np.zeros((6, 6))
        K_viga[1:3, 1:3] += self.EI * np.array([[12 / self.L ** 3, 6 / self.L ** 2],
                                                [6 / self.L ** 2, 4 / self.L]])
        K_viga[1:3, 4:6] += self.EI * np.array([[-12 / self.L ** 3, 6 / self.L ** 2],
                                                [-6 / self.L ** 2, 2 / self.L]])
        K_viga[4:6, 1:3] += self.EI * np.array([[-12 / self.L ** 3, -6 / self.L ** 2],
                                                [6 / self.L ** 2, 2 / self.L]])
        K_viga[4:6, 4:6] += self.EI * np.array([[12 / self.L ** 3, -6 / self.L ** 2],
                                                [-6 / self.L ** 2, 4 / self.L]])
        self.K = K_trelica + K_viga
        self.K = np.round(self.R @ self.K @ self.R.T, 2)
        self.K = sp.Matrix(self.K)

    def calcularMatrizdeForcas(self):
        """
        Compute element force vector.
        """
        x = sp.symbols('x')
        L = self.L
        p1, p2 = self.p1, self.p2
        q1, q2 = self.q1, self.q2

        # truss linear interpolation and stiffness matrix

        phi_a = (1-x/L)
        phi_b = x/L

        q = q1 + (q2-q1) * x / L

        f_trelica = sp.Matrix(np.zeros((6, 1)))
        f_trelica[0] += sp.integrate(phi_a * q, (x, 0, L))
        f_trelica[3] += sp.integrate(phi_b * q, (x, 0, L))

        f_viga = sp.Matrix(np.zeros((6, 1)))
        p = p1 + (p2 - p1) * x / L

        # Beam cubic interpolation and stiffness Matrix

        phi1 = (2 / L ** 3) * x ** 3 - (3 / L ** 2) * x ** 2 + 1
        phi2 = (1 / L ** 2) * x ** 3 - (2 / L) * x ** 2 + x
        phi3 = (-2 / L ** 3) * x ** 3 + (3 / L ** 2) * x ** 2
        phi4 = ((1 / L ** 2) * x ** 3 - (1 / L) * x ** 2)

        f_viga[1] += sp.integrate(phi1 * p, (x, 0, L))
        f_viga[2] += sp.integrate(phi2 * p, (x, 0, L))
        f_viga[4] += sp.integrate(phi3 * p, (x, 0, L))
        f_viga[5] += sp.integrate(phi4 * p, (x, 0, L))


        self.f = self.R * (f_viga + f_trelica)


