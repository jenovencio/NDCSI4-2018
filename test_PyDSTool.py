from PyDSTool import *
import matplotlib.pyplot as plt

"""
Multiple scales analysis of a duffing eq.
Excited near primary external resonance
See eq. 3.188 in 'Vibration and Stability', Jon J Thomsen
"""

# plot settings
plt.close('all')

# fixed simulation parameters
maxNumPoints = int(1e4)
maxStepSize = 1e-2
stepSize    = 1e-3
minStepSize = 1e-5

# fixed physical parameters
O = 0.0   # starting point for excitation
o = 1.0   # omega0 - linear eigenvalue
g = 0.5   # gamma - nonlinear parameter
q = 0.2   # q - external forcing level
b = 0.05  # beta - damping

# fixed equilibrium curve
# use that O = o + eps*s <=> s = O-o
ustr = '-b*o*u+q/(2*o)*sin(v)'  # amplitude
vstr = '(O-o)*u-3*g/(8*o)*u*u*u+q/(2*o)*cos(v)'  # phase

# predictor-corrector arguments
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['O']
PCargs.StepSize = stepSize
PCargs.MaxNumPoints = maxNumPoints
PCargs.MaxStepSize = maxStepSize
PCargs.MinStepSize = minStepSize
PCargs.LocBifPoints = ['LP', 'BP', 'B']
PCargs.StopAtPoints = ['B']
PCargs.verbosity = 2
PCargs.SaveEigen = True
# Note: 'MoorePenrose' is much more stable than'Natural' corrector:
PCargs.Corrector = 'MoorePenrose'

DSargs = args(name='duffing')
DSargs.pdomain = {'O': [0, 4.0]}  # domain boundaries
DSargs.pars = {'O': O, 'o': o, 'g': g, 'q': q, 'b': b}
DSargs.varspecs = {'u': ustr, 'v': vstr}
DSargs.ics = {'u': 0.01, 'v': 1.0}
testDS = Generator.Vode_ODEsystem(DSargs)

# loop beta's
bvec = [0.20, 0.10, 0.05]

fig1 = plt.figure(num=1)
for b in bvec:
    testDS.set(pars={'b': b})
    PyCont = ContClass(testDS)

    PyCont.newCurve(PCargs)

    print('Computing equilibrium curve...')
    start = clock()
    PyCont['EQ1'].forward()
    print('done in %.3f seconds!' % (clock()-start))
    PyCont.display(('O', 'u'), stability=True, linewidth=0.5, label='beta=%f' % b)