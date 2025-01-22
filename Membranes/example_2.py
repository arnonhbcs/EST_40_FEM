from Membranes.Node import No
from CST import CST
from Membranes.Plate import Placa

No1 = No(0, 0, 1, 'apoio_fixo')
No4 = No(0, 0.25, 4, 'apoio_fixo')
No2 = No(0.5, 0, 2,  'comum')
No3 = No(0.5, 0.25, 3,  'comum')

t = 0.025
E = 210e9
v = 0.3

elemento2 = CST([No1, No2, No3], t, E, v)
elemento1 = CST([No1, No3, No4], t, E, v)

elemento2.inserirPressaoEntreNos(2, 3, 7500e3)

nos = [No1, No2, No3, No4]
elementos = [elemento1, elemento2]

placa = Placa(nos, elementos)
placa.imporCondicoesContornos()
placa.resolverSistema()
placa.posProcessamento()

print('Vetor de Deslocamentos')
print(placa.d)

print('Tensões Principais de Cada Elemento')
for elemento in placa.elementos:
    print(elemento.calcularTensoesPrincipais())