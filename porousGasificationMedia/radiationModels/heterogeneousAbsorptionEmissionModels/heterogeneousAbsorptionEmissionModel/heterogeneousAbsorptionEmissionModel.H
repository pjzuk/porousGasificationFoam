/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiationModels::heterogeneousAbsorptionEmissionModel

Description
    Model to supply absorption and emission coefficients for
    heterogeneous radiation modelling

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousAbsorptionEmissionModel_H
#define heterogeneousAbsorptionEmissionModel_H

#include "absorptionEmissionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
             Class heterogeneousAbsorptionEmissionModel Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousAbsorptionEmissionModel 
:
    public Foam::radiation::absorptionEmissionModel
{

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("heterogeneousAbsorptionEmissionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heterogeneousAbsorptionEmissionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        heterogeneousAbsorptionEmissionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<heterogeneousAbsorptionEmissionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~heterogeneousAbsorptionEmissionModel();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient (net)
                virtual tmp<volScalarField> as(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> asCont(const label bandI = 0) const;

                //- Absorption coefficient (net)
                virtual tmp<volScalarField> borderAs(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                virtual tmp<volScalarField> borderAsCont(const label bandI = 0) const;


            // Emission coefficient

                //- emission coefficient (net)
                virtual tmp<volScalarField> es(const label bandi = 0) const;

                //- emission coefficient (net)
                virtual tmp<volScalarField> esCont(const label bandi = 0) const;

                //- emission coefficient (net)
                virtual tmp<volScalarField> borderEs(const label bandi = 0) const;

                //- emission coefficient (net)
                virtual tmp<volScalarField> borderEsCont(const label bandi = 0) const;


            // Border thickness layer
            virtual dimensionedScalar borderL(const label bandI = 0) const;

            //- Add radiative intensity for ray i
            virtual tmp<volScalarField> addIntensity
            (
                const label rayI,
                const volScalarField& ILambda
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
