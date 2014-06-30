/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::turbulentMapped2PatchFieldBase

Description
    Functionality for sampling fields using mappedPatchBase. Every call to
    mappedField() returns a sampled field, optionally scaled to maintain an
    area-weighted average.

    Example usage:

    {
        fieldName           T;          // default is same as fvPatchField
        setAverage          false;
        average             1.0;        // only if setAverage=true
        interpolationScheme cellPoint;  // default is cell
    }

SourceFiles
    turbulentMapped2PatchFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMapped2PatchFieldBase_H
#define turbulentMapped2PatchFieldBase_H

#include "fixedValueFvPatchFields.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mappedPatchBase;
template<class> class interpolation;

/*---------------------------------------------------------------------------*\
                  Class turbulentMapped2PatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class turbulentMapped2PatchFieldBase
{

protected:

    // Protected data

        //- Mapping engine
        const mappedPatchBase& mapper_;

        //- Underlying patch field
        const fvPatchField<Type>& patchField_;

        //- Name of field to sample
        word fieldName_;

        //- If true adjust the mapped field to maintain average value average_
        const bool setAverage_;

        //- Average value the mapped field is adjusted to maintain if
        //  setAverage_ is set true
        const Type average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;


public:

    // Constructors

        //- Construct from components
        turbulentMapped2PatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct from dictionary
        turbulentMapped2PatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const dictionary& dict
        );

        //- Construct empty
        turbulentMapped2PatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField
        );

        //- Construct copy
        turbulentMapped2PatchFieldBase
        (
            const turbulentMapped2PatchFieldBase<Type>& mapper
        );

        //- Construct copy, resetting patch and field
        turbulentMapped2PatchFieldBase
        (
            const mappedPatchBase& mapper,
            const fvPatchField<Type>& patchField,
            const turbulentMapped2PatchFieldBase<Type>& base
        );


    //- Destructor
    virtual ~turbulentMapped2PatchFieldBase<Type>()
    {}


    // Member functions

        //- Field to sample. Either on my or nbr mesh
        const GeometricField<Type, fvPatchField, volMesh>& sampleField() const;

        //- Map sampleField onto *this patch
        virtual tmp<Field<Type> > mappedField() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "turbulentMapped2PatchFieldBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
