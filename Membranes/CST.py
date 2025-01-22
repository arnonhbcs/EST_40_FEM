import numpy as np


class CST:
    def __init__(self, nos, t, E, v):
        """
        Implements a Constant Strain Triangle type of plate element.
        """
        self.matrizB = None
        self.D = None
        self.nos = nos
        self.A = nos[0]
        self.B = nos[1]
        self.C = nos[2]
        self.t = t
        self.E = E
        self.v = v
        self.K = None

        self.calcularK()
        self.F = np.zeros((6, 1))
        self.sigma = np.zeros((3, 1))

    def calcularK(self):
        """
        Computes stiffness matrix.
        """
        E = self.E
        v = self.v

        self.D = ((E / (1 - v ** 2)) *
                  np.array([
                      [1, v, 0],
                      [v, 1, 0],
                      [0, 0, (1 - v) / 2]
                  ]))

        A = self.A
        B = self.B
        C = self.C

        Area = (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y)) / 2

        self.matrizB = np.array([
            [B.y - C.y, 0, C.y - A.y, 0, A.y - B.y, 0],
            [0, C.x - B.x, 0, A.x - C.x, 0, B.x - A.x],
            [C.x - B.x, B.y - C.y, A.x - C.x, C.y - A.y, B.x - A.x, A.y - B.y]
        ]) / (2 * Area)

        self.K = self.matrizB.T @ self.D @ self.matrizB * self.t * Area

    def inserirPressaoEntreNos(self, x, y, p):
        """
        Insert pressure loads between nodes.
        """
        X = self.nos[x - 1]
        Y = self.nos[y - 1]
        L = np.sqrt((X.x - Y.x) ** 2 + (X.y - Y.y) ** 2)
        alpha = np.arctan2(Y.y - X.y, Y.x - X.x)
        s = np.sin(alpha)
        c = np.cos(alpha)

        self.F[2 * x - 2, 0] += p * L / 2 * s
        self.F[2 * x - 1, 0] += - p * L / 2 * c

        self.F[2 * y - 2, 0] += p * L / 2 * s
        self.F[2 * y - 1, 0] += - p * L / 2 * c
        self.F = np.round(self.F, 2)

    def calcularTensoesPrincipais(self):
        """
        Computes Von Mises stress.
        """
        sigmaMax = (self.sigma[0, 0] + self.sigma[1, 0]) / 2 + np.sqrt(
            ((self.sigma[0, 0] - self.sigma[1, 0]) / 2) ** 2 + self.sigma[2, 0] ** 2
        )
        sigmaMin = (self.sigma[0, 0] + self.sigma[1, 0]) / 2 - np.sqrt(
            ((self.sigma[0, 0] - self.sigma[1, 0]) / 2) ** 2 + self.sigma[2, 0] ** 2
        )

        return sigmaMax, sigmaMin
