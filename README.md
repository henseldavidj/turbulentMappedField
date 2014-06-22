turbulentMappedField
====================

re-implementation of mappedPatchFieldBase including scaling based on standard deviation of sampled field (OpenFOAM)

This is the early prototype of a mappedPatch boundary condition.  It has a very basic implementation at this time, but should work to scale a vector field's z-component based on a provided standard deviation and mean velocity.

I commented out the NoRepository business at the bottom of the header files, because this was causing the compiler to complain about redefinition
