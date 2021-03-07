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
    Foam::solidOde

Description
    An solid ODE solver for solid and gaseous chemistry

SourceFiles
    solidOde.C

\*---------------------------------------------------------------------------*/

#ifndef solidOde_H
#define solidOde_H

#include "solidChemistrySolver.H"
#include "chemistrySolver.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class SolidThermo, class GasThermo>
class solidOde;

/*---------------------------------------------------------------------------*\
                            Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class SolidThermo, class GasThermo>
class solidOde
:
    public solidChemistrySolver<CompType, SolidThermo, GasThermo>
{
    // Private data

        dictionary coeffsDict_;
        const word solverName_;
        autoPtr<ODESolver> odeSolver_;

        // Model constants

            scalar eps_;
            scalar scale_;


public:

    //- Runtime type information
    TypeName("solidOde");


    // Constructors

        //- Construct from components
        solidOde
        (
            ODESolidHeterogeneousChemistryModel<CompType, SolidThermo, GasThermo>& model,
            const word& modelName
        );


    //- Destructor
    virtual ~solidOde();


    // Member Functions

        scalar solve
        (
            scalarField& c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solidOde.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
