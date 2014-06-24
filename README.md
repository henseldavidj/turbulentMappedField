turbulentMappedField
====================

re-implementation of mappedPatchFieldBase including scaling based on standard deviation of sampled field (OpenFOAM)

This is the early prototype of a mappedPatch boundary condition.  It has a very basic implementation at this time, but should work to scale a vector field's z-component based on a provided (hard coded at the moment) standard deviation and mean velocity.  (See turbulentMappedPatchFieldBase.C for implementation)

Compiles on the Northeastern University Discovery Research Cluster (RHEL 6.3 Linux 64 bit)
