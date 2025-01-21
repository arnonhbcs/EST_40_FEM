import numpy as np
from Node import No
from Element import Elemento
from Frame import Portico

E = 200 * 10 ** 9 # Young Modulus
# Structure dimensions
d = 400 * 10 ** (-3)
b = 140 * 10 ** (-3)
tw = 6 * 10 ** (-3)
tf = 9 * 10 ** (-3)
h = d - 2 * tf

# area
A = h * tw + 2 * tf * b
# Stiffness
EA = E * A
I = 1.24204e-4
EI = E * I
# i
nos = []
elementos = []
num1 = 3  # numero de elementos finitos antes do apoio
num2 = 3  # numero de elementos finitos depois do apoio
x1 = np.linspace(0, 5.5, num1 + 1)
x2 = np.linspace(5.5, 5.5 + 2.5, num2 + 1)
x = np.concatenate((x1, x2[1:]))


# linear load distribution
def distribuicaoCarga(x):
    p1 = -90 * 10 ** 3
    p2 = -30 * 10 ** 3
    L = 5.5 + 2.5
    return p1 + (p2 - p1) * x / L


for i in range(len(x)):
    if x[i] == 0:
        tipo = 'engaste'
    elif x[i] == 5.5:
        tipo = 'apoio_livre_x'
    else:
        tipo = 'no_comum'
    nos.append(No(x[i], 0, i + 1, tipo))

for i in range(len(x) - 1):
    elementos.append(Elemento(A=nos[i], B=nos[i + 1], EA=EA, EI=EI, Id=i + 1,
                              p1=distribuicaoCarga(x[i]), p2=distribuicaoCarga(x[i + 1])))

portico = Portico(nos, elementos)
portico.calcularParametrosNodais()
portico.posProcessamento()
filename = 'exemplo_1.txt'
portico.escreverSolucao(filename)
