import numpy as np
import sympy as sp


class Placa:
    """
    Implements a discrete Membrane Structure with Finite Elements of the CST type (Constant Strain Triangle).
    """
    def __init__(self, nos, elementos):
        self.nos = nos
        self.elementos = elementos
        self.size = 2 * len(nos)
        self.K = sp.Matrix(np.zeros((self.size, self.size)))
        self.F = sp.Matrix(np.zeros((self.size, 1)))
        self.u = sp.Matrix(np.zeros((self.size, 1)))
        self.d = None

        self.simbolos = []
        self.calcularKeF()

    def calcularKeF(self):
        """
        Computes stiffness matrix and force vector
        """
        for elemento in self.elementos:
            A, B, C = elemento.A, elemento.B, elemento.C
            a, b, c = A.Id, B.Id, C.Id

            self.K[2 * a - 2: 2 * a, 2 * a - 2: 2 * a] += elemento.K[0:2, 0:2]
            self.K[2 * a - 2: 2 * a, 2 * b - 2: 2 * b] += elemento.K[0:2, 2:4]
            self.K[2 * a - 2: 2 * a, 2 * c - 2: 2 * c] += elemento.K[0:2, 4:6]

            self.K[2 * b - 2: 2 * b, 2 * a - 2: 2 * a] += elemento.K[2:4, 0:2]
            self.K[2 * b - 2: 2 * b, 2 * b - 2: 2 * b] += elemento.K[2:4, 2:4]
            self.K[2 * b - 2: 2 * b, 2 * c - 2: 2 * c] += elemento.K[2:4, 4:6]

            self.K[2 * c - 2: 2 * c, 2 * a - 2: 2 * a] += elemento.K[4:6, 0:2]
            self.K[2 * c - 2: 2 * c, 2 * b - 2: 2 * b] += elemento.K[4:6, 2:4]
            self.K[2 * c - 2: 2 * c, 2 * c - 2: 2 * c] += elemento.K[4:6, 4:6]

            self.F[2 * a - 2: 2 * a, 0] += elemento.F[0:2, 0].reshape(-1, 1)
            self.F[2 * b - 2: 2 * b, 0] += elemento.F[2:4, 0].reshape(-1, 1)
            self.F[2 * c - 2: 2 * c, 0] += elemento.F[4:6, 0].reshape(-1, 1)

    def inserirForcaConcentrada(self, Id, Fx, Fy):
        self.F[2 * Id - 2, 0] += Fx
        self.F[2 * Id - 1, 0] += Fy

    def imporCondicoesContornos(self):
        """
        Sets constraints.
        """
        for no in self.nos:
            Fxi = sp.symbols(f'Fx{no.Id}')
            Fyi = sp.symbols(f'Fy{no.Id}')
            uxi = sp.symbols(f'ux{no.Id}')
            uyi = sp.symbols(f'uy{no.Id}')

            if no.tipo == 'apoio_fixo':
                self.simbolos.append(Fxi)
                self.simbolos.append(Fyi)
                self.F[2 * no.Id - 2] += Fxi
                self.F[2 * no.Id - 1] += Fyi
            elif no.tipo == 'apoio_livre_x':
                self.simbolos.append(Fyi)
                self.simbolos.append(uxi)
                self.F[2 * no.Id - 1] += Fyi
                self.u[2 * no.Id - 2] += uxi
            elif no.tipo == 'apoio_livre_y':
                self.simbolos.append(Fxi)
                self.simbolos.append(uyi)
                self.F[2 * no.Id - 2] += Fxi
                self.u[2 * no.Id - 1] += uyi
            elif no.tipo == 'comum':
                self.simbolos.append(uyi)
                self.simbolos.append(uxi)
                self.u[2 * no.Id - 2] += uxi
                self.u[2 * no.Id - 1] += uyi

    def resolverSistema(self):
        """
        Compute node displacements
        """
        sol = sp.solve(self.K * self.u - self.F, self.simbolos)
        disp = {str(key): sol[key] for key in sol.keys() if str(key).startswith('u')}
        d = np.zeros((self.size, 1))

        for key in disp.keys():
            Id = int(''.join([char for char in str(key) if char.isdigit()]))
            if str(key).startswith('ux'):
                d[2 * Id - 2, 0] = disp[key]
            elif str(key).startswith('uy'):
                d[2 * Id - 1, 0] = disp[key]
        self.d = d

    def posProcessamento(self):
        for elemento in self.elementos:
            a, b, c = elemento.A.Id, elemento.B.Id, elemento.C.Id
            d = np.zeros((6, 1))
            d[0, 0] = self.d[2 * a - 2, 0]
            d[1, 0] = self.d[2 * a - 1, 0]
            d[2, 0] = self.d[2 * b - 2, 0]
            d[3, 0] = self.d[2 * b - 1, 0]
            d[4, 0] = self.d[2 * c - 2, 0]
            d[5, 0] = self.d[2 * c - 1, 0]

            elemento.sigma = elemento.D @ elemento.matrizB @ d
