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

\*---------------------------------------------------------------------------*/

#include "turbulentMappedPatchFieldBase.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
turbulentMappedPatchFieldBase<Type>::turbulentMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const word& fieldName,
    const bool setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average),
    interpolationScheme_(interpolationScheme)
{}


template<class Type>
turbulentMappedPatchFieldBase<Type>::turbulentMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_
    (
        dict.template lookupOrDefault<word>
        (
            "fieldName",
            patchField_.dimensionedInternalField().name()
        )
    ),
    setAverage_(readBool(dict.lookup("setAverage"))),
    average_(pTraits<Type>(dict.lookup("average"))),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.lookup("interpolationScheme") >> interpolationScheme_;
    }
}


template<class Type>
turbulentMappedPatchFieldBase<Type>::turbulentMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(patchField_.dimensionedInternalField().name()),
    setAverage_(false),
    average_(pTraits<Type>::zero),
    interpolationScheme_(interpolationCell<Type>::typeName)
{}


template<class Type>
turbulentMappedPatchFieldBase<Type>::turbulentMappedPatchFieldBase
(
    const turbulentMappedPatchFieldBase<Type>& mapper
)
:
    mapper_(mapper.mapper_),
    patchField_(mapper.patchField_),
    fieldName_(mapper.fieldName_),
    setAverage_(mapper.setAverage_),
    average_(mapper.average_),
    interpolationScheme_(mapper.interpolationScheme_)
{}


template<class Type>
turbulentMappedPatchFieldBase<Type>::turbulentMappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const turbulentMappedPatchFieldBase<Type>& base
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(base.fieldName_),
    setAverage_(base.setAverage_),
    average_(base.average_),
    interpolationScheme_(base.interpolationScheme_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const GeometricField<Type, fvPatchField, volMesh>&
turbulentMappedPatchFieldBase<Type>::sampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    if (mapper_.sameRegion())
    {
        if (fieldName_ == patchField_.dimensionedInternalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    patchField_.dimensionedInternalField()
                );
        }
        else
        {
            const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
            return thisMesh.template lookupObject<fieldType>(fieldName_);
        }
    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(fieldName_);
    }
}


template<class Type>
tmp<Field<Type> > turbulentMappedPatchFieldBase<Type>::mappedField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());

    // Result of obtaining remote values
    tmp<Field<Type> > tnewValues(new Field<Type>(0));
    Field<Type>& newValues = tnewValues();

    switch (mapper_.mode())
    {
        case mappedPatchBase::NEARESTCELL:
        {
            const mapDistribute& distMap = mapper_.map();

            if (interpolationScheme_ != interpolationCell<Type>::typeName)
            {
                // Send back sample points to the processor that holds the cell
                vectorField samples(mapper_.samplePoints());
                distMap.reverseDistribute
                (
                    (
                        mapper_.sameRegion()
                      ? thisMesh.nCells()
                      : nbrMesh.nCells()
                    ),
                    point::max,
                    samples
                );

                autoPtr<interpolation<Type> > interpolator
                (
                    interpolation<Type>::New
                    (
                        interpolationScheme_,
                        sampleField()
                    )
                );
                const interpolation<Type>& interp = interpolator();

                newValues.setSize(samples.size(), pTraits<Type>::max);
                forAll(samples, cellI)
                {
                    if (samples[cellI] != point::max)
                    {
                        newValues[cellI] = interp.interpolate
                        (
                            samples[cellI],
                            cellI
                        );
                    }
                }
            }
            else
            {
                newValues = sampleField();
            }

            distMap.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        case mappedPatchBase::NEARESTPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorIn
                (
                    "void turbulentMappedPatchFieldBase<Type>::updateCoeffs()"
                )<< "Unable to find sample patch " << mapper_.samplePatch()
                 << " in region " << mapper_.sampleRegion()
                 << " for patch " << patchField_.patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            mapper_.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), pTraits<Type>::zero);

            const fieldType& nbrField = sampleField();

            forAll(nbrField.boundaryField(), patchI)
            {
                const fvPatchField<Type>& pf =
                    nbrField.boundaryField()[patchI];
                label faceStart = pf.patch().start();

                forAll(pf, faceI)
                {
                    allValues[faceStart++] = pf[faceI];
                }
            }

            mapper_.distribute(allValues);
            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "turbulentMappedPatchFieldBase<Type>::updateCoeffs()"
            )<< "Unknown sampling mode: " << mapper_.mode()
             << nl << abort(FatalError);
        }
    }

    if (setAverage_)
    {
/*
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());
*/
            //Info<< "averagePsi: "<< averagePsi << endl;
	
	// Take streamwise component of sampled field "newValues"
	// Changed to X component (0), as typical streamwise coordinate
        scalarField valZ(newValues.component(0));
        //Info<< "valZ: "<< valZ << endl;
   
        //TODO: parameterize this:
        // Desired mean velocity.  This should be provided at the boundary condition dictionary level
        scalar meanVelocityWanted = 61;

        //TODO: Parameterize this
        // Desired standard deviation.  This should be provided at the boundary condition dictionary level
        
        // Area weighted fluctuation, from Sydney data. As percentage of mean velocity 61. (7.72307) based on
        // third-order polynomial fit to data series, and calculated by reimann sum like 'integration' 12.66% 
        scalar stdevGiven = 0.1266 * meanVelocityWanted;

        // this scalar field might not be necessary, but didn't like the average_ parameter.
        scalarField avgGiven(newValues.size(), meanVelocityWanted);

	// Calculate the normalized field by dividing the sampled field by the field's average
        scalarField normZ = valZ / avgGiven;
        //Info<< "normZ "<< normZ <<endl;

	// NOTE: As derived on 7/24/14, it is not necessary to scale a normalized field, rather
	// any field will do.  You just need the standard deviation of the sampled field to create
	// a scaling factor.
	
	// Calculate the average of the normalized field
        scalar avgZNorm =
                    gSum(patchField_.patch().magSf()*normZ)
                   /gSum(patchField_.patch().magSf());
        //Info<< "avgZNorm "<< avgZNorm <<endl;

	// Calculate the standard deviation of the normalized field
	scalarField diffFromMean = avgZNorm - normZ;
	//Info<< "diffFromMean: "<< diffFromMean << endl;
	scalarField diffsqr = sqr(diffFromMean);
	//Info<< "diffsqr: "<< diffsqr << endl;
	scalar var = gSum(patchField_.patch().magSf()*diffsqr)/gSum(patchField_.patch().magSf());
	//Info<< "var: "<< var << endl;
	scalar stdev = sqrt(var);
	//Info<< "stdev "<< stdev << endl;

	// Avoid divide by zero error
	scalar sdScaleFactor = 1;
	if (stdev != 0){
	sdScaleFactor = stdevGiven / stdev;
	}
	//Info << "sdScaleFactor " << sdScaleFactor << endl;

	// Scale normalized field by multipling by the standard deviation scaling field
	// This provides a field with the desired standard deviation
	scalarField sdScaled = normZ * sdScaleFactor;
	normZ *= sdScaleFactor;
	//Info<< "sdScaled "<< normZ << endl;

	// Calculate mean value of newly scaled field
	scalar sdscaledMean =
	                    gSum(patchField_.patch().magSf()*normZ)
	                   /gSum(patchField_.patch().magSf());

	// Calculate the difference between the desired mean and the newly scaled field mean value
	scalar meanShift = sdscaledMean - meanVelocityWanted;

	// Shift the newly scaled values by the difference calculated previously ("meanShift")
	sdScaled -= meanShift;
	//Info<< "scaled mean shifted values"<< sdScaled << endl;

	// Reassign this field back to boundary by replacing components of newValues
	// Cross-stream components are assumed zero, based on initial conditions
	
	//newValues.component(2) = shiftedValues;
	newValues.replace(0,sdScaled);  //changed scaled component to X component
	newValues.replace(1,0);
	newValues.replace(2,0);
	//Info<< "scaled newValues "<< newValues << endl;


		// unit test:
		scalar meanValueTest = gSum(patchField_.patch().magSf()*newValues.component(0))
				                   /gSum(patchField_.patch().magSf());
		Info<< "test meanValue: "<< meanValueTest << endl;
	
		scalarField diffFromMean2 = meanValueTest - newValues.component(0);
	
		scalarField diffsqr2 = sqr(diffFromMean2);
				//Info<< "diffsqr: "<< diffsqr << endl;
	
		scalar var2 = gSum(patchField_.patch().magSf()*diffsqr2)/gSum(patchField_.patch().magSf());
		//Info<< "var: "<< var << endl;
	
		scalar stdev2 = sqrt(var2);
		Info<< "stdev given: " << stdevGiven << " = stdev of test: " << stdev2 << endl;


/*
//Not needed for my scaling

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
*/
    }


    // Restore tag
    UPstream::msgType() = oldTag;

//writes out entire mapped field
    //Info<< "newValues: "<< newValues << endl;
    return tnewValues;
}


template<class Type>
void turbulentMappedPatchFieldBase<Type>::write(Ostream& os) const
{
    os.writeKeyword("fieldName") << fieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
    os.writeKeyword("interpolationScheme") << interpolationScheme_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
