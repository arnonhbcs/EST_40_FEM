class No:
    """
    Structure Node
    """

    def __init__(self, x, y, Id, tipo='comum', x_prescrito=0, y_prescrito=0):
        # Node coordinates
        self.x = x
        self.y = y
        # Node key
        self.Id = Id
        # Node displacements
        self.ux = 0
        self.vy = 0
        self.theta = 0
        # type of constraint (translation, rotation)
        self.tipo = tipo
        # offset
        self.x_prescrito = x_prescrito
        self.y_prescrito = y_prescrito
        # Node forces
        self.Fx = 0
        self.Fy = 0
        self.M = 0


