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

#include "coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
//#include "regionProperties.H"
#include "basicThermo.H"
#include "LESModel.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::
coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    neighbourFieldName_("undefined-neighbourFieldName"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    KName_("undefined-K"),
    QrIncident_(0.0),
    distanceFromBase_(0.0),
    angleOfIncidence_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::
coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
(
    const coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    radiationCoupledBase
    (
        p,
        ptf.emissivityMethod(),
        ptf.emissivity_
    ),
    neighbourFieldName_(ptf.neighbourFieldName_),
    neighbourFieldRadiativeName_(ptf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(ptf.fieldRadiativeName_),
    KName_(ptf.KName_),
    QrIncident_(ptf.QrIncident_),
    distanceFromBase_(ptf.distanceFromBase_),
    angleOfIncidence_(ptf.angleOfIncidence_)
//    emissivity_(ptf.emissivity_)
{}


coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::
coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    radiationCoupledBase(p, dict),
    neighbourFieldName_(dict.lookup("neighbourFieldName")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    KName_(dict.lookup("K")),
    QrIncident_(readScalar(dict.lookup("QrIncident"))),
    distanceFromBase_(readScalar(dict.lookup("distanceFromBase"))),
    angleOfIncidence_(readScalar(dict.lookup("angleOfIncidence"))),
    viewFactors_(viewFactorFunction())
{
    /*if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::"
            "coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
*/
    //Info << "viewFactors: " << viewFactors_ << endl;
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::
coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
(
    const coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField&
        wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    radiationCoupledBase
    (
        wtcsf.patch(),
        wtcsf.emissivityMethod(),
        wtcsf.emissivity_
    ),
    neighbourFieldName_(wtcsf.neighbourFieldName_),
    neighbourFieldRadiativeName_(wtcsf.neighbourFieldRadiativeName_),
    fieldRadiativeName_(wtcsf.fieldRadiativeName_),
    KName_(wtcsf.KName_),
    QrIncident_(wtcsf.QrIncident_),
    distanceFromBase_(wtcsf.distanceFromBase_),
    angleOfIncidence_(wtcsf.angleOfIncidence_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::K() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

/*    if (KName_ == "none")
    {
        const compressible::LESModel& model =
        db().lookupObject<compressible::LESModel>("LESProperties");
        
        tmp<volScalarField> talpha = model.alphaEff();

        const basicThermo& thermo =
            db().lookupObject<basicThermo>("thermophysicalProperties");

        return
            talpha().boundaryField()[patch().index()]
           *thermo.Cp()().boundaryField()[patch().index()];
    }*/

    if (KName_ == "solidThermo")
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>("thermophysicalProperties");

        //scalarField K_ = thermo.kappa(patch().index());
        return thermo.kappa(patch().index());
    }

    else if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
    {
        return patch().lookupPatchField<volScalarField, scalar>(KName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KName_))
    {
        const symmTensorField& KWall =
            patch().lookupPatchField<volSymmTensorField, scalar>(KName_);

        vectorField n(patch().nf());

        return n & KWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField"
            "::K() const"
        )   << "Did not find field " << KName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'K' to 'none', 'solidThermo', a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    //const mapDistribute& distMap = mpp.map();
    //const mappedPatchBase& distMap = mpp;

    scalarField intFld(patchInternalField());

    const coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField&
    nbrField =  refCast
        <const coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarList radField(nbrPatch.size(),0.0);  
    scalarField Twall(patch().size(),0.0);

    // In solid
    if(neighbourFieldRadiativeName_ != "none") //nbr Radiation Qr
    {
        radField =
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldRadiativeName_
            );
   
        // Swap to obtain full local values of neighbour radiative heat flux field
        mpp.distribute(radField);

//        emissivity_ =
//            patch().lookupPatchField<volScalarField, scalar>("emissivity");

//        const scalarField temissivity = emissivity();
            const fvMesh& mesh = patch().boundaryMesh().mesh();
            const radiation::radiationModel& radiation =
                mesh.lookupObject<radiation::radiationModel>
                (
                    "radiationProperties"
                );

            scalarField temissivity
            (
                radiation.absorptionEmission().e()().boundaryField()
                [
                    //nbrPatch.index()
                    patch().index()
                ]
            );

        scalarField nbrTotalFlux(-radField);

/*
        Twall =
            (radField + myKDelta*intFld + nbrKDelta*nbrIntFld)
            /(myKDelta + nbrKDelta);
*/

        //hard code for now.
        const scalar sigma_ = 5.67e-8;

        scalarField gradT_((((QrIncident_*viewFactors_)+radField)*temissivity-temissivity*sigma_*pow(intFld,4))/K());

        forAll(*this, i)
        {
		//fixed gradient BC, use internal T to replace Tb for re-radiation
                this->refValue()[i] = operator[](i);  // not used
                this->refGrad()[i] = gradT_[i];
                this->valueFraction()[i] = 0.0;
        }

    }
    else // In fluid
    {
        radField =
            patch().lookupPatchField<volScalarField, scalar>
            (
                fieldRadiativeName_
            );
	this->operator==(this->patchInternalField());
//        Twall = nbrIntFld;

  /*      this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0; */
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qr = gSum(radField*patch().magSf());
	if (neighbourFieldRadiativeName_ != "None"){ Info << "Solid: ";}
	else {Info << "Fluid: ";}

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            //<< this->dimensionedInternalField().name() << " -> "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            //<< this->dimensionedInternalField().name() << " :"
            << " radiativeFlux:" << Qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}

    
scalarField coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::viewFactorFunction() const
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
            scalar angle1 = Foam::acos(zz/sqrt(ghost2faceDistances[i][j]));    // hard coded to positivez direction
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
    


void coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< fieldRadiativeName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("QrIncident")<< QrIncident_
        << token::END_STATEMENT << nl;
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivityMode") << emissivityMethod() << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
