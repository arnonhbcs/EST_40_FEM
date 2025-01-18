class No:
    """
    Essa classe implementa um nó de elemento finito.
    """

    def __init__(self, x, y, Id, tipo='comum', x_prescrito=0, y_prescrito=0):
        self.x = x
        self.y = y
        self.Id = Id  # número do nó
        self.ux = 0
        self.vy = 0
        self.theta = 0
        self.tipo = tipo
        self.x_prescrito = x_prescrito
        self.y_prescrito = y_prescrito
        self.Fx = 0
        self.Fy = 0
        self.M = 0


