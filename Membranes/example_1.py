from Node import No
from CST import CST
from Plate import Placa

No1 = No(x=0, y=0, Id=1, tipo='apoio_fixo')
No2 = No(x=3, y=0, Id=2, tipo='apoio_livre_x')
No3 = No(x=3/2, y=2.5/2, Id=3, tipo='comum')
No4 = No(x=3, y=2.5, Id=4, tipo='comum')
No5 = No(x=0, y=2.5, Id=5, tipo='comum')
nos = [No1, No2, No3, No4, No5]

elemento2 = CST(nos=[No2, No4, No3], t=8e-3, E=105e9, v=0.35)
# Verify if the given load is N/m or pressure.
t = 8e-3
E = 105e9
v = 0.35
elemento2.inserirPressaoEntreNos(2, 1, 15e6 * t)

elementos = [
    CST(nos=[No1, No2, No3], t=t, E=E, v=v), elemento2,
    CST(nos=[No4, No5, No3], t=t, E=E, v=v), CST(nos=[No5, No1, No3], t=t, E=E, v=v),
]
placa = Placa(nos=nos, elementos=elementos)
placa.inserirForcaConcentrada(5, 10e3, -20e3)
placa.imporCondicoesContornos()
placa.resolverSistema()
placa.posProcessamento()


