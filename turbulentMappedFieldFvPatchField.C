/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "turbulentMappedFieldFvPatchField.H"

#include "volFields.H"
#include "interpolationCell.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    turbulentMappedPatchFieldBase<Type>(*this, *this)
{}


template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const turbulentMappedFieldFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    turbulentMappedPatchFieldBase<Type>(*this, *this, ptf)
{}


template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    turbulentMappedPatchFieldBase<Type>(*this, *this, dict)
{}


template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,

    // mappedPatchBase
    const word& sampleRegion,
    const sampleMode sampleMode,
    const word& samplePatch,
    const scalar distance,

    // My settings
    const word& fieldName,
    const bool setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase
    (
        p.patch(),
        sampleRegion,
        sampleMode,
        samplePatch,
        distance
    ),
    turbulentMappedPatchFieldBase<Type>
    (
        *this,
        *this,
        fieldName,
        setAverage,
        average,
        interpolationScheme
    )
{}


template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const turbulentMappedFieldFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    turbulentMappedPatchFieldBase<Type>(ptf)
{}


template<class Type>
turbulentMappedFieldFvPatchField<Type>::turbulentMappedFieldFvPatchField
(
    const turbulentMappedFieldFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    turbulentMappedPatchFieldBase<Type>(*this, *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void turbulentMappedFieldFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

	this->operator==(this->turbulentMappedPatchFieldBase<Type>::mappedField());

    if (debug)
    {
        Info<< "operating on field:" << this->dimensionedInternalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void turbulentMappedFieldFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    turbulentMappedPatchFieldBase<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
