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
	
	// Take components of sampled field
    	scalarField sampleX(newValues.component(0));
	scalarField sampleY(newValues.component(1));
	scalarField sampleZ(newValues.component(2));
	//Info<< "sampleX: "<< sampleX << endl;
	//Info<< "sampleY: "<< sampleY << endl;
    	//Info<< "sampleZ: "<< sampleZ << endl;
   
    	//TODO: parameterize this:
   	// Desired mean velocity.  This should be provided at the boundary condition dictionary level
   	// Changed to 45.32604 per area-weighted average velocity from curve fit of data
   	// Changed to centerline velocity, as this corresponds to what's shown in the second two results sheets.
   	scalar meanVelocityWantedX = 70.76;
	scalar meanVelocityWantedY = 0;
	scalar meanVelocityWantedZ = 0;

    	//TODO: Parameterize this
    	// Desired standard deviation.  This should be provided at the boundary condition dictionary level
        
    	// Area weighted fluctuation, from Sydney data. As percentage of mean velocity 61. (7.80) based on
    	// fourth-order polynomial fit to fit data series, and calculated by reimann sum like 'integration' 12.66% 
    
	// modification to average of 1% fluctuations, in u',v',w'

    	scalar stdevGivenX = 0.01 * meanVelocityWantedX;
	scalar stdevGivenY = 0.01 * meanVelocityWantedX; // Assumed v' and w' are same as u'
	scalar stdevGivenZ = 0.01 * meanVelocityWantedX; // Assumed v' and w' are same as u'

    // this scalar field might not be necessary, but didn't like the average_ parameter.
    scalarField avgGiven(newValues.size(), meanVelocityWantedX);

	// Calculate the average of the normalized field
        scalar sampleXMean =
                    gSum(patchField_.patch().magSf()*sampleX)
                   /gSum(patchField_.patch().magSf());

		scalar sampleYMean =
					gSum(patchField_.patch().magSf()*sampleY)
					/gSum(patchField_.patch().magSf());

		scalar sampleZMean =
					gSum(patchField_.patch().magSf()*sampleZ)
					/gSum(patchField_.patch().magSf());
        //Info<< "sampleXMean: "<< sampleXMean <<endl;
		//Info<< "sampleYMean: "<< sampleYMean <<endl;
		//Info<< "sampleZMean: "<< sampleZMean <<endl;

	// Calculate the standard deviation of the sampled X component
	scalarField diffFromMeanX = sampleXMean - sampleX;
	//Info<< "diffFromMean: "<< diffFromMean << endl;
	scalarField diffsqrX = sqr(diffFromMeanX);
	//Info<< "diffsqr: "<< diffsqr << endl;
	scalar varX = gSum(patchField_.patch().magSf()*diffsqrX)/gSum(patchField_.patch().magSf());
	//Info<< "var: "<< var << endl;
	scalar stdevX = sqrt(varX);
	//Info<< "stdev "<< stdev << endl;

	// Avoid divide by zero error
	scalar sdScaleFactorX = 1;
	if (stdevX != 0){
	sdScaleFactorX = stdevGivenX / stdevX;
	}
	//Info << "sdScaleFactor " << sdScaleFactor << endl;

	// Calculate the standard deviation of the sampled Y component
	scalarField diffFromMeanY = sampleYMean - sampleY;
	scalarField diffsqrY = sqr(diffFromMeanY);
	scalar varY = gSum(patchField_.patch().magSf()*diffsqrY)/gSum(patchField_.patch().magSf());
	scalar stdevY = sqrt(varY);

	// Avoid divide by zero error
	scalar sdScaleFactorY = 1;
	if (stdevY != 0){
	sdScaleFactorY = stdevGivenY / stdevY;
	}
	//Info << "sdScaleFactor " << sdScaleFactor << endl;

	// Calculate the standard deviation of the sampled Z component
	scalarField diffFromMeanZ = sampleZMean - sampleZ;
	scalarField diffsqrZ = sqr(diffFromMeanZ);
	scalar varZ = gSum(patchField_.patch().magSf()*diffsqrZ)/gSum(patchField_.patch().magSf());
	scalar stdevZ = sqrt(varZ);

	// Avoid divide by zero error
	scalar sdScaleFactorZ = 1;
	if (stdevZ != 0){
	sdScaleFactorZ = stdevGivenZ / stdevZ;
	}



	// Scale normalized field by multipling by the standard deviation scaling field
	// This provides a field with the desired standard deviation
	scalarField sdScaledX = sampleX * sdScaleFactorX;
	sampleX *= sdScaleFactorX;
	//Info<< "sdScaled "<< normX << endl;

	// Calculate mean value of newly scaled field
	scalar sdscaledMeanX =
	                    gSum(patchField_.patch().magSf()*sampleX)
	                   /gSum(patchField_.patch().magSf());

	// Calculate the difference between the desired mean and the newly scaled field mean value
	scalar meanShiftX = sdscaledMeanX - meanVelocityWantedX;

	// Shift the newly scaled values by the difference calculated previously ("meanShift")
	sdScaledX -= meanShiftX;
	//Info<< "scaled mean shifted values"<< sdScaled << endl;

	// Scale sampleY field
	scalarField sdScaledY = sampleY * sdScaleFactorY;
	sampleY *= sdScaleFactorY;

	// Calculate mean value of newly scaled field
	scalar sdscaledMeanY =
	                    gSum(patchField_.patch().magSf()*sampleY)
	                   /gSum(patchField_.patch().magSf());

	// Calculate the difference between the desired mean and the newly scaled field mean value
	scalar meanShiftY = sdscaledMeanY - meanVelocityWantedY;
	// Shift the newly scaled values by the difference calculated previously ("meanShift")
	sdScaledY -= meanShiftY;
	//Info<< "scaled mean shifted values"<< sdScaled << endl;

	// Scale sampleZ field
	scalarField sdScaledZ = sampleZ * sdScaleFactorZ;
	sampleZ *= sdScaleFactorZ;
	//Info<< "sdScaled "<< normX << endl;

	// Calculate mean value of newly scaled field
	scalar sdscaledMeanZ =
	                    gSum(patchField_.patch().magSf()*sampleZ)
	                   /gSum(patchField_.patch().magSf());

	// Calculate the difference between the desired mean and the newly scaled field mean value
	scalar meanShiftZ = sdscaledMeanZ - meanVelocityWantedZ;
	// Shift the newly scaled values by the difference calculated previously ("meanShift")
	sdScaledZ -= meanShiftZ;
	//Info<< "scaled mean shifted values"<< sdScaled << endl;
	

	// Reassign this field back to boundary by replacing components of newValues
	// Cross-stream components are assumed zero, based on initial conditions
	
	//newValues.component(2) = shiftedValues;
	newValues.replace(0,sdScaledX);  //changed scaled component to X component
	newValues.replace(1,sdScaledY);
	newValues.replace(2,sdScaledZ);
	//Info<< "scaled newValues "<< newValues << endl;


		// unit test:
		scalar meanValueTestX = gSum(patchField_.patch().magSf()*newValues.component(0))
				                   /gSum(patchField_.patch().magSf());
		Info<< "Scaled mean X "<< meanValueTestX << endl;
	
		scalarField diffFromMeanTestX = meanValueTestX - newValues.component(0);
	
		scalarField diffsqrTestX = sqr(diffFromMeanTestX);
				//Info<< "diffsqr: "<< diffsqr << endl;
	
		scalar varTestX = gSum(patchField_.patch().magSf()*diffsqrTestX)/gSum(patchField_.patch().magSf());
		//Info<< "var: "<< var << endl;
	
		scalar stdevTestX = sqrt(varTestX);
		Info<< "Scaled stDev X: " << stdevTestX << endl;

		// unit test:
		scalar meanValueTestY = gSum(patchField_.patch().magSf()*newValues.component(1))
				                   /gSum(patchField_.patch().magSf());
		Info<< "Scaled mean Y "<< meanValueTestY << endl;
	
		scalarField diffFromMeanTestY = meanValueTestY - newValues.component(1);
	
		scalarField diffsqrTestY = sqr(diffFromMeanTestY);
				//Info<< "diffsqr: "<< diffsqr << endl;
	
		scalar varTestY = gSum(patchField_.patch().magSf()*diffsqrTestY)/gSum(patchField_.patch().magSf());
		//Info<< "var: "<< var << endl;
	
		scalar stdevTestY = sqrt(varTestY);
		Info<< "Scaled stDev Y: " << stdevTestY << endl;

		// unit test:
		scalar meanValueTestZ = gSum(patchField_.patch().magSf()*newValues.component(2))
				                   /gSum(patchField_.patch().magSf());
		Info<< "Scaled mean Z "<< meanValueTestZ << endl;
	
		scalarField diffFromMeanTestZ = meanValueTestZ - newValues.component(2);
	
		scalarField diffsqrTestZ = sqr(diffFromMeanTestZ);
				//Info<< "diffsqr: "<< diffsqr << endl;
	
		scalar varTestZ = gSum(patchField_.patch().magSf()*diffsqrTestZ)/gSum(patchField_.patch().magSf());
		//Info<< "var: "<< var << endl;
	
		scalar stdevTestZ = sqrt(varTestZ);
		Info<< "Scaled stDev Z: " << stdevTestZ << endl;

		
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
