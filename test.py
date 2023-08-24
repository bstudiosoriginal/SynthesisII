from initialconditions import Input, IntermediateAngles, Output, PositionDisplacement
from system import System

R1 = complex(0, 0)
R2 = complex(-1.4, -0.76)
R3 = complex(-1, -2.3)

# initial conditions
displacement = PositionDisplacement(R1, R2, R3)
input = Input(PHI1=0, PHI2=126, PHI3=252)
intermediate = IntermediateAngles(gam1=0, gam2=-6, gam3=37)
output = Output(PSI1=0, PSI2=33, PSI3=37)

# system solver
sys = System(displacements=displacement, _input=input, _intermediate=intermediate, _output=output)
sys.solve()

# show
sys.show()

# plot
sys.plot()
