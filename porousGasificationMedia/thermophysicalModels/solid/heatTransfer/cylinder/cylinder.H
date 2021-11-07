/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

Description
    Cylinder heat transfer model

SourceFiles
    cylinder.C

\*---------------------------------------------------------------------------*/

#ifndef cylinder_H
#define cylinder_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

class cylinderCONV
:
    public heatTransferModel
{
private:

    scalar cylinderRadius_;
    const volVectorField& Up_;
    const volScalarField& rhop_;
    const volScalarField& alphap_;
    const volScalarField& mup_;
    const basicThermo& thermop_;

public:

    //- Runtime type information
    TypeName("cylinderCONV");


    // Constructors

        //- Construct from components
        cylinderCONV
        (
            const volScalarField& por,
	        const volScalarField& por0
        );

    // Selectors

        //- Return a reference to the selected heatTransfer model
        static autoPtr<cylinderCONV> New
        (
            const volScalarField& por,
	        const volScalarField& por0
        );


    //- Destructor
    virtual ~cylinderCONV()
    {}


    // Member Functions

        //- Return the field for heat transfer between solid/gas on surfaces
        virtual tmp<volScalarField> CONV() const;

        //- Read heatTransferProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
