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

#include "myFlowRateFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myFlowRateFvPatchScalarField::
myFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("none"),
    massFluxFraction_(1.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::myFlowRateFvPatchScalarField::            // This constructor is called
myFlowRateFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    massFluxFraction_(dict.lookupOrDefault<scalar>("massFluxFraction", 1.0))
{
    //Info<< "debug1:This Constructor called. "<< endl;

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

Foam::myFlowRateFvPatchScalarField::
myFlowRateFvPatchScalarField
(
    const myFlowRateFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    massFluxFraction_(ptf.massFluxFraction_)
{}


Foam::myFlowRateFvPatchScalarField::
myFlowRateFvPatchScalarField
(
    const myFlowRateFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}

Foam::myFlowRateFvPatchScalarField::
myFlowRateFvPatchScalarField
(
    const myFlowRateFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myFlowRateFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::myFlowRateFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


void Foam::myFlowRateFvPatchScalarField::updateCoeffs()
{
    /*Info<< "debug1:Update coeffs() run. "<< endl;
    Info << "This before: "<< *this <<endl;
    Info << "Value Fraction before: "<< valueFraction()<<endl;
    Info << "refvalue before: "<< refValue() <<endl;*/
    if (this->updated())
    {
        return;
    }

    const label patchi = patch().index();

    const LESModel<EddyDiffusivity<compressible::turbulenceModel>>& turbModel =
        db().lookupObject
        <
            LESModel<EddyDiffusivity<compressible::turbulenceModel>>
        >
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);         // this is the value of phi (kg/s) during the time step at the patch

    const scalarField alphap(turbModel.alphaEff(patchi));       // effective thermal diffusivity for enthalpy for a patch (kg/m.s)
    //Info << "phi: "<<phip<< endl;
    refValue() = massFluxFraction_;
    refGrad() = 0.0;
    
    valueFraction() =                                   // update value fraction at the patch
        1.0
        /
        (
            1.0 +
            alphap*patch().deltaCoeffs()*patch().magSf()/max(mag(phip), SMALL)
        );          // deltaCoeffs retuns face2cell distance
                    //  magSf <surface area at patch of each mesh cell>
                    // mag(phip) returns absolute values of the vectors
    mixedFvPatchField<scalar>::updateCoeffs();
    /*Info<< "This after! "<< *this<< endl;
    Info << "Value Fraction after : "<< valueFraction()<<endl;
    Info << "refvalue after: "<< refValue() <<endl;*/
    if (debug)
    {
        scalar phi = gSum(-phip*(*this));
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }
}


void Foam::myFlowRateFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
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
        myFlowRateFvPatchScalarField
    );

}

// ************************************************************************* //
