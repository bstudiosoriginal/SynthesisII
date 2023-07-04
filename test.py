from initialconditions import Input, IntermediateAngles, Output, PositionDisplacement
from system import System

R1 = complex(0.1, 0.1)
R2 = complex(-0.07+0.1, 0.4+0.1)
R3 = complex(-0.3+0.1, 0.7+0.1)

# initial conditions
displacement = PositionDisplacement(R1, R2, R3)
input = Input(PHI1=0, PHI2=50, PHI3=75)
intermediate = IntermediateAngles(gam1=0, gam2=7, gam3=12)
output = Output(PSI1=0, PSI2=22.5, PSI3=45)

# system solver
sys = System(displacements=displacement, _input=input, _intermediate=intermediate, _output=output)
sys.solve()

# show
sys.show()

# plot
sys.plot()
