import sympy


class Radius():

    def __init__(self, R, theta) -> None:
        self.Rx = R * sympy.cos(theta)
        self.Ry = R * sympy.sin(theta)

class MassObject():

    def __init__(self, mass, R, L) -> None:
        self.mass = mass
        self.R = R
        self.L = L


class System():
    
    def __init__(self, *mass_objects):
        self.mass_objects = mass_objects

    def solve(self):
        total = 0
        for i in range(len(self.mass_objects)):
            mass = self.mass_objects[i].mass
            r = self.mass_objects[i].R
            l = self.mass_objects[i].L
            solution0 = mass * r.Rx * l
            solution1 = mass * r.Ry * l


