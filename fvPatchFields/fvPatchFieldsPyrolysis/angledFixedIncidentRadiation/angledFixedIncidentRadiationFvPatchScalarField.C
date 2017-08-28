/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "angledFixedIncidentRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "constants.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "RectangularMatrix.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

angledFixedIncidentRadiationFvPatchScalarField::
angledFixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    QrIncident_(p.size(), 0.0),
    distanceFromBase_(0.0),
    angleOfIncidence_(0.0)
{}


angledFixedIncidentRadiationFvPatchScalarField::
angledFixedIncidentRadiationFvPatchScalarField
(
    const angledFixedIncidentRadiationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    QrIncident_(psf.QrIncident_),
    distanceFromBase_(psf.distanceFromBase_),
    angleOfIncidence_(psf.angleOfIncidence_)
{}


angledFixedIncidentRadiationFvPatchScalarField::
angledFixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    QrIncident_("QrIncident", dict, p.size()),
    distanceFromBase_(readScalar(dict.lookup("distanceFromBase"))),
    angleOfIncidence_(readScalar(dict.lookup("angleOfIncidence"))),
    viewFactors_(viewFactorFunction())
{
    //Info << "viewFactors: " << viewFactors_ << endl;
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


    
    
angledFixedIncidentRadiationFvPatchScalarField::
angledFixedIncidentRadiationFvPatchScalarField
(
    const angledFixedIncidentRadiationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    QrIncident_(psf.QrIncident_),
    distanceFromBase_(psf.distanceFromBase_),
    angleOfIncidence_(psf.angleOfIncidence_)
{}


angledFixedIncidentRadiationFvPatchScalarField::
angledFixedIncidentRadiationFvPatchScalarField
(
    const angledFixedIncidentRadiationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    temperatureCoupledBase(patch(), ptf),
    QrIncident_(ptf.QrIncident_),
    distanceFromBase_(ptf.distanceFromBase_),
    angleOfIncidence_(ptf.angleOfIncidence_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void angledFixedIncidentRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    QrIncident_.autoMap(m);
}


void angledFixedIncidentRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(psf, addr);

    const angledFixedIncidentRadiationFvPatchScalarField& thftpsf =
        refCast<const angledFixedIncidentRadiationFvPatchScalarField>
        (
            psf
        );

    QrIncident_.rmap(thftpsf.QrIncident_, addr);
}


void angledFixedIncidentRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    
    scalarField intFld(patchInternalField());
    
    const radiation::radiationModel& radiation =
        db().lookupObject<radiation::radiationModel>("radiationProperties");

            scalarField temissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );

    // gradient()[i] to access, intFld[i] to access
    gradient() =
        temissivity*
        (
            QrIncident_*viewFactors_ - physicoChemical::sigma.value()*pow4(intFld)  //use internal cell T
        )/kappa(*this);
    
    
    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qr = gSum(kappa(*this)*gradient()*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " -> "
            << " radiativeFlux:" << Qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}

scalarField angledFixedIncidentRadiationFvPatchScalarField::viewFactorFunction() const
    {
    // - COMPUTATION OF VIEW FACTORS at angles - //
    // - TODO: include distance feature so that fluxer can be located at a height

    // compute position from where flux is incident
    const vectorField base(patch().Cf());
    vectorField locationOfFluxer(base);
    // apply rotation matrix to obtain final position of ghost-patch
    scalar angle = constant::mathematical::pi*angleOfIncidence_/180.0;
    scalar cosine = Foam::cos(angle);
    scalar sine = Foam::sin(angle);
    tensor rotMat(1, 0, 0, 0, cosine, -sine, 0, sine, cosine);
    forAll (patch(), i){
        locationOfFluxer[i] = rotMat & locationOfFluxer[i];
    }

    scalar m = locationOfFluxer.size();
    scalar n = patch().size();
        
    // Now exists cell centers for the patch and the ghost-incident radiator cell centers
    RectangularMatrix<scalar> ghost2faceDistances(m, n);    // store distance between fluxer-faces and patch-faces
    RectangularMatrix<scalar> beta1(m, n);  // store patch viewfactor angles
    RectangularMatrix<scalar> beta2(m, n);  // store fluxer viewfactor angles
    
    forAll (patch(), j){
        forAll (locationOfFluxer, i){       // assume each face center as frame of reference to compute distances
            scalar xx = locationOfFluxer[i].x() - base[j].x();
            scalar yy = locationOfFluxer[i].y() - base[j].y();
            scalar zz = locationOfFluxer[i].z() - base[j].z();
            // Unable to reference to face normals to compute angles. Hence hard coded below to compute angle1
            // ideally i want to be able to
            // multiply by surface normals so that angle is computed automatically
            ghost2faceDistances[i][j] = sqr(xx) + sqr(yy) + sqr(zz);
            // Vector<scalar> v(xx,yy,zz);
            // initialize beta values as zeros
            beta1[j][i] = 0;
            beta2[i][j] = 0;
            //if (mag(xx) <= SMALL){
            scalar angle1 = Foam::acos(zz/sqrt(ghost2faceDistances[i][j]));    // hard coded to positive z direction
                                                                                // will not work if angle provided is 0
                    
            beta2[i][j] = Foam::cos(angle - angle1)*patch().magSf()[i];
            beta1[j][i] = Foam::cos(angle1)*patch().magSf()[i]
                                            /
                    (constant::mathematical::pi*ghost2faceDistances[i][j]);
           //     }
        }
    }
   // Info << "beta1" <<beta1 << endl;
   // Info << "beta2" <<beta2 << endl;
    RectangularMatrix<scalar> shape(m, n);
    shape = (1/patch().magSf()[0])*beta1*beta2;
    //Info << "shape" << shape << endl;
    scalarField viewFactor(patchInternalField());
    forAll(patch(), j){
        forAll(locationOfFluxer,i){
            if (i==j){
                if (shape[i][j] > 1){ shape[i][j] = 1;}
                viewFactor[i] = shape[i][j];
            }
        }
    }

    return viewFactor;
}

void angledFixedIncidentRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    temperatureCoupledBase::write(os);
    QrIncident_.writeEntry("QrIncident", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    angledFixedIncidentRadiationFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
