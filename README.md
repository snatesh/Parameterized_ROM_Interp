# Parameterized_ROM_Interp
Local reduced order basis interpolation applied to the parabolic diffusion equation with a random coefficient field

In this work, a local reduced order basis (RB) interpolation method is described and explored numerically on a parameterized parabolic diffusion equation with a log-normal coefficient field. RBs and their corresponding reduced order models (ROMs) are often not robust with respect to changes in the parameter λ ∈ R^p to which the RB, generated by proper orthogonal decomposition (POD) of transient solution approximations, is associated. As such, interpolation of a set of RBs corresponding to a collection of parameters Λ = {λ_i} may be a sensible approach to approximating the RB for parameter λ not in Λ. However, straightforward interpolation of RBs does not guarantee that the interpolated vectors satisfy the properties of a RB, namely linear independence and orthogonality. The interpolation method considered here, and first introduced by Amsallem and Farhat, represents a RB generated by POD and associated to a parameter, or “operating point”, as a point on the Grassmanian manifold. Representing RBs in this way offers a more amenable starting point to their interpolation, and results from differential geometry
actually expose an algorithm to do so. 

This repository constains prototype MATLAB codes that apply the local RB interpolation alogirthm on the parabolic PDE and assess its efficacy in terms of a standard cross-validation procedure. We use rectangular finite elements and a tensor product construction to discretize the problem in 2D, and the Karhuenen-Loeve expansion is used to efficiently sample the log-normal random field on the grid.

In the PDF document, MethodDescription.pdf, the considered model problem is presented, its discretization reviewed, and the local RB interpolation algorithm is explicated and applied on the model problem. Lastly, an alternate formulation that overcomes the computational bottleneck inherent to the local RB interpolation method for problems with non-affine parameter dependence is described.

The following animation illustrates how the parabolic diffusion equation with a random coefficient field "randomly" and rapidly smooths initial noisy data during the transient phase constrained by the boundary conditions, towards the equilibrium distribution which is depicted in Figure 4 of the summary pdf.

![Sample Frame](transient_smooth.gif)

