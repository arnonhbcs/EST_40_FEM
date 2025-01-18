import numpy as np
import sympy as sp
from Element import Elemento
from Node import No
from typing import List


class Portico:
    """
    Implementation of a frame structure for calculating forces
    and displacements.
    """
    def __init__(self, nos: List[No], elementos: List[Elemento]):
        self.nos = nos
        self.size = 3 * len(self.nos)
        self.elementos = elementos
        self.K = np.zeros((self.size, self.size))
        self.f = np.zeros((self.size, 1))
        self.solution = None
        self.esforcos = None
        self.calcularMatrizRigidez()
        self.calcularMatrizdeForcas()

    def calcularMatrizRigidez(self):
        for elemento in self.elementos:
            a, b = elemento.A.Id, elemento.B.Id
            self.K[3 * a - 3:3 * a, 3 * a - 3:3 * a] += elemento.K[0:3, 0:3]
            self.K[3 * a - 3:3 * a, 3 * b - 3:3 * b] += elemento.K[0:3, 3:6]
            self.K[3 * b - 3:3 * b, 3 * a - 3:3 * a] += elemento.K[3:6, 0:3]
            self.K[3 * b - 3:3 * b, 3 * b - 3:3 * b] += elemento.K[3:6, 3:6]

    def calcularMatrizdeForcas(self):
        for elemento in self.elementos:
            a, b = elemento.A.Id, elemento.B.Id
            self.f[3 * a - 3:3 * a, 0:1] += elemento.f[0:3, 0:1].reshape(3, 1)
            self.f[3 * b - 3:3 * b, 0:1] += elemento.f[3:6, 0:1].reshape(3, 1)

    def inserirEsforcosExternos(self, F, M, Id):
        """
        Insere esforços externos em nó.
        :param Id: número do nó.
        :param F: Vetor de Força Concentrada
        :param M: Momento em um nó.
        """
        self.f[3 * Id - 3, 0] += F[0]
        self.f[3 * Id - 2, 0] += F[1]
        self.f[3 * Id - 1, 0] += M

    def calcularParametrosNodais(self):
        """
        Essa função resolve o sistema Ku = f para determinar deslocamentos
        e giros
        :return: dict com a solução
        """
        K_cc = self.K
        # impor condicoes de contorno
        for no in self.nos:
            linha_zerada = np.zeros((1, self.size))
            coluna_zerada = np.zeros((self.size, 1))
            Id = no.Id
            if no.tipo == 'engaste':
                K_cc[3 * Id - 3, :] = linha_zerada
                K_cc[:, 3 * Id - 3] = coluna_zerada.flatten()
                K_cc[3 * Id - 3, 3 * Id - 3] = 1
                self.f[3 * Id - 3, 0] = 0

                K_cc[3 * Id - 2, :] = linha_zerada
                K_cc[:, 3 * Id - 2] = coluna_zerada.flatten()
                K_cc[3 * Id - 2, 3 * Id - 2] = 1
                self.f[3 * Id - 2, 0] = 0

                K_cc[3 * Id - 1, :] = linha_zerada
                K_cc[:, 3 * Id - 1] = coluna_zerada.flatten()
                K_cc[3 * Id - 1, 3 * Id - 1] = 1
                self.f[3 * Id - 1, 0] = 0

            elif no.tipo == 'apoio_livre_x':
                K_cc[3 * Id - 2, :] = linha_zerada
                K_cc[:, 3 * Id - 2] = coluna_zerada.flatten()
                K_cc[3 * Id - 2, 3 * Id - 2] = 1
                self.f[3 * Id - 2, 0] = no.y_prescrito

            elif no.tipo == 'apoio_livre_y':
                K_cc[3 * Id - 3, :] = linha_zerada
                K_cc[:, 3 * Id - 3] = coluna_zerada.flatten()
                K_cc[3 * Id - 3, 3 * Id - 3] = 1
                self.f[3 * Id - 3, 0] = 0

            elif no.tipo == 'apoio_duplo_y_prescrito':

                K_cc[3 * Id - 3, :] = linha_zerada
                K_cc[:, 3 * Id - 3] = coluna_zerada.flatten()
                K_cc[3 * Id - 3, 3 * Id - 3] = 1
                self.f[3 * Id - 3, 0] = 0

                K_cc[3 * Id - 2, :] = linha_zerada
                K_cc[3 * Id - 2, 3 * Id - 2] = 1
                self.f[3 * Id - 2, 0] = no.y_prescrito

            elif no.tipo == 'apoio_duplo_x_prescrito':
                K_cc[3 * Id - 2, :] = linha_zerada
                K_cc[:, 3 * Id - 2] = coluna_zerada.flatten()
                K_cc[3 * Id - 2, 3 * Id - 2] = 1
                self.f[3 * Id - 2, 0] = 0

                K_cc[3 * Id - 3, :] = linha_zerada
                K_cc[3 * Id - 3, 3 * Id - 3] = 1
                self.f[3 * Id - 3, 0] = no.x_prescrito

        # resolver sistema
        self.solution = np.linalg.inv(K_cc) @ self.f

    def posProcessamento(self):
        u = self.solution
        u = np.array(u, dtype='float64')[:, np.newaxis]

        for elemento in self.elementos:
            a, b = elemento.A.Id, elemento.B.Id

            uE = np.vstack((u[3 * a - 3:3 * a, 0], u[3 * b - 3:3 * b, 0]))
            uL = elemento.R.T @ uE

            x = sp.symbols('x')
            L = elemento.L
            phi1 = (2 / L ** 3) * x ** 3 - (3 / L ** 2) * x ** 2 + 1
            phi2 = (1 / L ** 2) * x ** 3 - (2 / L) * x ** 2 + x
            phi3 = (-2 / L ** 3) * x ** 3 + (3 / L ** 2) * x ** 2
            phi4 = ((1 / L ** 2) * x ** 3 - (1 / L) * x ** 2)

            v = uL[1, 0] * phi1 + uL[2, 0] * phi2 + uL[4, 0] * phi3 + uL[5, 0] * phi4

            elemento.N = elemento.EA * (uL[3, 0] - uL[0, 0]) / L
            elemento.V = -elemento.EI * v.diff(x, 3)
            elemento.M = elemento.EI * v.diff(x, 2).subs({x: L / 2})

    def escreverSolucao(self, filename):
        parametrosNodais = self.solution
        with open(filename, 'w') as f:
            f.write('Parametros Nodais: \n \n')
            for i in range(len(parametrosNodais)):
                if i % 3 == 0:
                    f.write(f'\n u{i + 1} = {parametrosNodais[i, 0]} m \n')
                elif i % 3 == 1:
                    f.write(f' v{i + 1} = {parametrosNodais[i, 0]} m \n')
                elif i % 3 == 2:
                    f.write(f' theta{i + 1} = {parametrosNodais[i, 0]} rad \n')

            for i in range(len(self.elementos)):
                f.write(f'\n Elemento {i + 1} \n')
                f.write(f'N = {round(self.elementos[i].N / 1000, 3)} kN \n')
                f.write(f'V = {round(self.elementos[i].V / 1000, 3)} kN \n')
                f.write(f'M(L/2) = {round(self.elementos[i].M / 1000, 3)} kN.m \n')
