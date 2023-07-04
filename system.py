import rich.table
import rich.console
import numpy as np


def det(m):
    # determinant of 2x2 matrix
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]

def plot_sol(V, origin, lim_x=[-20, 20], lim_y=[-20, 20], scale=1):
    import matplotlib.pyplot as plt
    plt.quiver(*origin, V[:, 0], V[:, 1], angles='xy', scale_units='xy', scale=scale)
    # plt.plot(V) 
    plt.xlim(lim_x)
    plt.ylim(lim_y)
    plt.show()

# constant
i = complex(0, 1)

class System():

    """
    Class that solves a four bar linkage problem
    """

    def __init__(self, displacements, _input, _intermediate, _output):
        self._displacements = displacements
        self._input = _input
        self._intermediate = _intermediate
        self._output = _output

    def solve(self):
        R1, R2, R3 = self._displacements.R1, self._displacements.R2, self._displacements.R3
        d2 = R2 - R1
        d3 = R3 - R1 
        # print(d2, d3)
        # change in gamma in radians
        ga1, ga2, ga3 = self._intermediate.gam1, self._intermediate.gam2, self._intermediate.gam3
        g2 = np.radians(ga2 - ga1) # which angle
        g3 = np.radians(ga3 - ga1) # angle diff of coupler link


        # LHS of equation
        phi1, phi2, phi3 = self._input.phi1, self._input.phi2, self._input.phi3
        ph2 = np.radians(phi2 - phi1)
        ph3 = np.radians(phi3 - phi1)

        # RHS of equation
        psi1, psi2, psi3 = self._output.psi1, self._output.psi2, self._output.psi3
        ps2 = np.radians(psi2 - psi1)
        ps3 = np.radians(psi3 - psi1)

        # A matrix for the system of equations
        A1 = np.array([[np.exp(ph2*i)-1, np.exp(g2*i)-1], [np.exp(ph3*i)-1, np.exp(g3*i)-1]])
        Bz2 = np.array([[d2, np.exp(g2*i)-1], [d3, np.exp(g3*i)-1]])
        Bz5 = np.array([[np.exp(ph2*i)-1, d2], [np.exp(ph3*i)-1, d3]])

        A2 = np.array([[np.exp(ps2*i)-1, np.exp(g2*i)-1], [np.exp(ps3*i)-1, np.exp(g3*i)-1]])
        Bz4 = np.array([[d2, np.exp(g2*i)-1], [d3, np.exp(g3*i)-1]])
        Bz6 = np.array([[np.exp(ps2*i)-1, d2], [np.exp(ps3*i)-1, d3]])

        Z4 = det(Bz4)/det(A2)
        Z6 = det(Bz6)/det(A2)
        Z5 = det(Bz5)/det(A1)
        Z2 = det(Bz2)/det(A1)
        Z3 = Z5-Z6
        Z1 = Z2 + Z3 - Z4

        Z = np.array([Z1, Z2, Z3, Z4, Z5, Z6])

        L = np.abs(Z) 
        th = np.degrees(np.angle(Z))
        X = np.real(Z)
        Y = np.imag(Z)
        rows = []
        for a, b, c, d, e in zip(Z.tolist(), X.tolist(), Y.tolist(), L.tolist(), th.tolist()):
            rows.append([a, b, c, d, e])
        table = rich.table.Table(title='4 bar Synthesis')
        columns = ['Z vector ', 'X position', 'Y position', 'Lenght', 'Angle (deg)']
        for column in columns:
            table.add_column(column)
        # print('Z', )
        for row in rows:
            nr = []
            for r in row:
                if isinstance(r, complex):
                    nr.append(str(round(r.real, 4))+' + j'+str(round(r.imag, 4)))
                else:
                    nr.append(str(round(r, 4)))
                # print(nr)
            table.add_row(*nr)
        self.table = table
        self.X = X
        self.Y = Y
        self.Z = Z
        self.L = L
        self.th = th

    def show(self):
        try:
            console = rich.console.Console()
            console.print(self.table)
        except Exception:
            print('You must solve first')
    
    def plot(self):
        # origins of the vectors
        try:
            X = self.X
            Y = self.Y
            Z1 = self.Z[0]
            Z2 = self.Z[1]
            Z3 = self.Z[2]
            Z4 = self.Z[3]
            Z5 = self.Z[4]
            Z6 = self.Z[5]
            l = np.array([[0, 0], [0, 0], [Z2.real, Z2.imag], [Z1.real, Z1.imag], [Z2.real, Z2.imag], [Z1.real+Z4.real, Z1.imag+Z4.imag]])

            # biggest point 
            MAX_POINT_X, MAX_POINT_Y = np.max([*(X+l.T[0]), 0]), np.max([*(Y+l.T[1]), 0])
            MIN_POINT_X, MIN_POINT_Y = np.min([*(X+l.T[0]), 0]), np.min([*(Y+l.T[1]), 0])

            MAX_POINT = np.max([MAX_POINT_X, MAX_POINT_Y])
            MIN_POINT = np.max([MIN_POINT_X, MIN_POINT_Y])
            # shift the origins
            # l -= np.array(np.abs(MAX_POINT_X-MIN_POINT_X)/2, np.abs(MAX_POINT_Y-MIN_POINT_Y)/2)
            # plot
            plot_sol(np.vstack([X, Y]).T, l.T, [MIN_POINT_X-5, MAX_POINT_X+5], [MIN_POINT_Y-5, MAX_POINT_Y+5])
        except Exception:
            print('You must solve first!!')
        
