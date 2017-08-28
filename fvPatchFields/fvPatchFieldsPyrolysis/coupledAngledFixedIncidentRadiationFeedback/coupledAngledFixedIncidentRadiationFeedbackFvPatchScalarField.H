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
    coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField

Description
    Mixed boundary condition for convection/radiation coupling to be used with
    chtMultiRegionRadiation solvers.

    Evaluate convective and radiative heat fluxes on each side.

    In the solid side if the heat flux is outgoing fixed value is used and if
    the heat flux is incoming fixedGradient is used.

    In the fluid side an average temperature is calculated and used as fixed
    value BC.

    Example usage:

    In the solid:
    upperWall_to_channel
    {
        type            compressible::coupledAngledFixedIncidentRadiationFeedback;
        value           uniform 500;
        neighbourFieldName T;
        neighbourFieldRadiativeName Qr;
        fieldRadiativeName none;
        K               K;
    }

    In the fluid:
    channel_to_upperWall
    {
        type            compressible::coupledAngledFixedIncidentRadiationFeedback;
        value           uniform 500;
        neighbourFieldName T;
        neighbourFieldRadiativeName none;
        fieldRadiativeName Qr;
        K               K;
    }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: if K is "none" looks up RASModel and basicThermo, otherwise expects
    the solver to calculate a 'K' field.

SourceFiles
    coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField_H
#define coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "radiationCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
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
    
        // - Distance from Base
        const scalar distanceFromBase_;
    
        // - Angle of incidence
        const scalar angleOfIncidence_;
    
        scalarField viewFactors_;

public:

    //- Runtime type information
    TypeName("compressible::coupledAngledFixedIncidentRadiationFeedback");


    // Constructors

        //- Construct from patch and internal field
        coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const
            coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
        (
            const coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField&,
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
                new coupledAngledFixedIncidentRadiationFeedbackFvPatchScalarField
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
    
        // return view factors
        scalarField viewFactorFunction() const;
    
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
