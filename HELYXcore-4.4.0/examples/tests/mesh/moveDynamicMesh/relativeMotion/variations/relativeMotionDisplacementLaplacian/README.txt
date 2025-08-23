Test case demonstrating the combination of a solid-body motion solver
with a Laplacian mesh distortion solver.

There is an outer cylinder, spinning with a constant rpm. Inside it,
a smaller cylinder counter-spins (with constant but opposite rpm).
Finally, inside this smaller cylinder there is a relative motion
between two blocks (using boundary conditions and a mesh distortion solver).
