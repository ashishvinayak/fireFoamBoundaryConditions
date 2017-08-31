/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::
    compressible::
    coupledFixedIncidentRadiationFeedbackFvPatchScalarField

Description
    Mixed boundary condition for convection/radiation coupling to be used with
    chtMultiRegionRadiation solvers.

    Evaluate convective and radiative heat fluxes on each side.

    In the solid side if the heat flux is outgoing fixed value is used and if
    the heat flux is incoming fixedGradient is used.

    In the fluid side an average temperature is calculated and used as fixed
    value BC.

    Radiant flux from fire is radiated back at the surface of material. 
 
    Example usage:

    In the solid:
    region0_to_pmmaRegion_pmma
    {
        type            compressible::coupledFixedIncidentRadiationFeedback;
        value           $internalField;
        neighbourFieldName            T;		// name of the field — no or T
        neighbourFieldRadiativeName	Qr;
        fieldRadiativeName		none;
        K		solidThermo;
        QrIncident	50000;                      // W/m2
        emissivityMode	solidRadiation;
    }
    In the fluid:
    coupledPatch
    {
        type            compressible::coupledFixedIncidentRadiationFeedback;
        value           $internalField;
        neighbourFieldName            T;		// name of the field — no or T
        neighbourFieldRadiativeName	none;
        fieldRadiativeName		Qr;
        K		K;
        emissivityMode	solidRadiation;
        QrIncident	0;
    }


    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: if K is "none" looks up RASModel and basicThermo, otherwise expects
    the solver to calculate a 'K' field.

SourceFiles
    coupledFixedIncidentRadiationFeedbackFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFixedIncidentRadiationFeedbackFvPatchScalarField_H
#define coupledFixedIncidentRadiationFeedbackFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "radiationCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class coupledFixedIncidentRadiationFeedbackFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class coupledFixedIncidentRadiationFeedbackFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public radiationCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

         //- Name of the radiative heat flux on the neighbout region
        const word neighbourFieldRadiativeName_;

        //- Name of the radiative heat flux on the my region
        const word fieldRadiativeName_;

        //- Name of thermal conductivity field
        const word KName_;

        //- Incident radiative heat flux
        const scalar QrIncident_;

        //- Emissivity
//        scalarField emissivity_;

public:

    //- Runtime type information
    TypeName("compressible::coupledFixedIncidentRadiationFeedback");


    // Constructors

        //- Construct from patch and internal field
        coupledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        coupledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const
            coupledFixedIncidentRadiationFeedbackFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledFixedIncidentRadiationFeedbackFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        coupledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const coupledFixedIncidentRadiationFeedbackFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledFixedIncidentRadiationFeedbackFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        tmp<scalarField> K() const;

        //- Get corresponding Qr field
        const scalarField& Qr() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
