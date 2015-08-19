# SPH_misc

Some code I wrote for SPH simulations. Very bare bone. 
The code quality/functionality is _very_ likely to improve in the near future, but expecting much before beginning of September would be a mistake. Before then look at your own risk and don't expect much!

Warning: few good code writing practices were followed during the development of these projects.

SPH_2D
- simple 2D SPH implementation 
- based on Alg 1 in STAR paper on SPH (Eurographics 2014)
https://graphics.ethz.ch/~sobarbar/papers/Sol14/2014_EG_SPH_STAR.pdf

PBF_2D 
- simple 2D implementation of PBF liquid simulation
- based on the paper by Macklin
http://mmacklin.com/pbf_sig_preprint.pdf

The following 3 projects are based on the paper by Takahashi:
http://gamma.cs.unc.edu/ViscousSPH/

Explicit_visc_2D
- 2D PBF with explicit full form viscosity 

Implicit_visc_2D
- 2D PBF with implicit full form viscosity 

Implicit_visc_3D 
- 3D PBF with implicit full form viscosity 
- likely something wrong with coefficients, since it doesn't look good/accurate
- Note: the world is rotating around Y axis in this one

Any questions/comments/concerns/feedback can be addressed to apazylbekova@gmail.com
