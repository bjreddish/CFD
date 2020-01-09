# 2D Supersonic and Hypersonic Viscous Flow Over a Flat Plate

This example looks at a viscous supersonic flow over a flat plate using the MacCormack's technique. This example is provided by chapter 10 of John D. Anderson's book Computational Fluid Dynamics: The Basics with Applications.  The flow passes over a sharp flat plate with no angle of attack. The plate of length L will form a laminar boundary layer at the beginning of the plate and the flow was assumed to remain laminar for the remainder of the length. Due to this boundary layer, the flow senses a curvature that generates a curved shock wave starting at the leading edge. The space between the plate and the shock wave is called the shock layer. This problem may be simple, but the results can provide insight into the complex phenomenon that occur in the shock layer along with the boundary layer.

Three versions of the code exist to model different conditions of flow. The first one, superVisc.py, is equipped to simulate supersonic flows. The python script hypervisc.py can simulate hypersonic flow over a flat plate. Lastly, the script interact.py simulates two parallel plates, one on the top of the domain and the other on the bottom. 

The scripts produce data for the flow field at steady state at the end of the simulation and save them as .h5 files. These files can then be processed with the post processing python scripts in the same directory. 

A brief overview of the theory used in the code, and some results are recorded in a [jupyter notebook]( https://nbviewer.jupyter.org/github/bjreddish/CFD/blob/master/superVisc/superViscFlatPlate.ipynb)


![Hypersonic Velocity](https://github.com/bjreddish/CFD/blob/master/superVisc/images/hypersonicVel.png)
