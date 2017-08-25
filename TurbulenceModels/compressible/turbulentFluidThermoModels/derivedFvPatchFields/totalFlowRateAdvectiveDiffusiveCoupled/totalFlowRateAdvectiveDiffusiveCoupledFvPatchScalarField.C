/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "turbulenceModel.H"
// Added by me
#include "mappedPatchBase.H"
#include "mapDistribute.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    nbrPhiName_("nbrPhi"),
    rhoName_("none"),
    massFluxFraction_(1.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}

// OH HELLO
Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::            // This constructor is called
totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    nbrPhiName_(dict.lookupOrDefault<word>("nbrPhi", "phiGas")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    massFluxFraction_(dict.lookupOrDefault<scalar>("massFluxFraction", 1.0))
{
    Info<< "debug1:This Constructor called. "<< endl;
    
    refValue() = 1.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
    
    if (dict.found("value"))
    {
        //Info<< "debug1:Value found. "<< endl;   // value found. print value
        fvPatchField<scalar>::operator=
        (
         Field<scalar>("value", dict, p.size())
         );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}

Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    rhoName_(ptf.rhoName_),
    massFluxFraction_(ptf.massFluxFraction_)
{}


Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    nbrPhiName_(tppsf.nbrPhiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}

Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    nbrPhiName_(tppsf.nbrPhiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


void Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::updateCoeffs()
{
    //Info << "updateCoeffs outer!" <<endl;
    if (this->updated())
    {
        //Info<< "Updated. "<< endl<< *this<<endl;
        return;
    }
    
    // get coupling information from mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];
    // look for phiGas in the neighbor region
    scalarList phin =
    -nbrPatch.lookupPatchField<surfaceScalarField, scalar>(nbrPhiName_);
    mpp.distribute(phin);
    
    const label patchi = patch().index();       // at the patch look for value of alphap

    const compressible::turbulenceModel& turbulence = db().lookupObject<compressible::turbulenceModel>("turbulenceModel");
    const scalarField alphap(turbulence.alphaEff(patchi));       // effective thermal diffusivity for enthalpy for a patch (kg/m.s)
    
    refValue() = 1;                             // compute value fraction acc to new_val = valueFrac * refValue
    refGrad() = 0.0;
    //Info <<"phin:"<<phin << endl;
    // ideally the iteration should look at each coupledpatchid to include multiple couplings
    forAll(*this, i){
        if (-phin[i] > 0) {
            this->valueFraction()[i] =  1.0;	//(1.0 + alphap[i]*patch().deltaCoeffs()[i]*patch().magSf()[i]/max(mag(phin[i]),SMALL));    // with eps usually
        }
            else { this->valueFraction()[i] = 0; }
    }

    mixedFvPatchField<scalar>::updateCoeffs();

 
    if (debug)
    {
        scalar phi = gSum(-phin*(*this));
       
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
    //        << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }
}


void Foam::totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("nbrPhi") << nbrPhiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("massFluxFraction") << massFluxFraction_
        << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalFlowRateAdvectiveDiffusiveCoupledFvPatchScalarField
    );

}

// ************************************************************************* //
