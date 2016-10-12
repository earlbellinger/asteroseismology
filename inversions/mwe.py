import adipy
import matplotlib.pyplot as pl

css, eigs = adipy.read_amde('modelS/modelS.amde')
fgong = adipy.load_fgong('modelS/modelS.fgong')
Ks = {'rho': adipy.kernels('rho', fgong, css, eigs), 
      'c': adipy.kernels('c', fgong, css, eigs)}
r = fgong['var'][::-1,0] # FGONG's are written from surface to centre; adipy knows this
R = fgong['glob'][1]
x = r/R

# make a plot similar to Gough (1991), Fig. 9

I = [i for i,cs in enumerate(css) if cs[17]==1 and cs[18]==21][0]
pl.plot(x, R*Ks['c'][I])
pl.plot(x, R*Ks['rho'][I])
pl.axis([-0.05,1.05,-5.5,5.5])
pl.show()
