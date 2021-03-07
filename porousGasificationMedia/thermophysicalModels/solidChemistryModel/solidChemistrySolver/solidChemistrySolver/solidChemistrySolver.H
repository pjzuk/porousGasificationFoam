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
    Foam::soliChemistrySolver

Description
    An abstract base class for solving chemistry

SourceFiles
    solidChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidChemistrySolver_H
#define solidChemistrySolver_H

#include "ODEChemistryModel.H"
#include "ODESolidHeterogeneousChemistryModel.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class SThermoType, class GThermoType>
class solidChemistrySolver;

/*---------------------------------------------------------------------------*\
                      Class solidChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class SThermoType, class GThermoType>
class solidChemistrySolver
{
protected:

    // Protected data

        //- Reference to the chemistry model
        ODESolidHeterogeneousChemistryModel<CompType, SThermoType, GThermoType>& model_;

        //- Name of the chemistry solver
        const word name_;


public:

        //- Runtime type information
        TypeName("solidChemistrySolver");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            solidChemistrySolver,
            dictionary,
            (
                ODESolidHeterogeneousChemistryModel<CompType, SThermoType, GThermoType>& model,
                const word& modelName
            ),
            (model, modelName)
        );


    // Constructors

        //- Construct from components
        solidChemistrySolver
        (
            ODESolidHeterogeneousChemistryModel<CompType, SThermoType, GThermoType>& model,
            const word& modelName
        );


    //- Selector
    static autoPtr<solidChemistrySolver> New
    (
        ODESolidHeterogeneousChemistryModel<CompType, SThermoType, GThermoType>& model,
        const word& compTypeName,
        const word& thermoTypeName
    );


    //- Destructor
    virtual ~solidChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolidChemistrySolver(Comp, SThermo, GThermo)                      \
                                                                              \
    typedef solidChemistrySolver<Comp, SThermo, GThermo>                      \
        solidChemistrySolver##Comp##SThermo##GThermo;                         \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        solidChemistrySolver##Comp##SThermo##GThermo,                         \
        "solidChemistryModel<"#Comp","#SThermo","#GThermo">",                 \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        solidChemistrySolver##Comp##SThermo##GThermo,                         \
        dictionary                                                            \
    );


#define makeSolidChemistrySolverType(SS,Comp, SThermo, GThermo)               \
                                                                              \
    typedef SS<Comp, SThermo, GThermo>                                        \
        SS##Comp##SThermo##GThermo;                                           \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##SThermo##GThermo,                                           \
        #SS"<"#Comp","#SThermo","#GThermo">",                                 \
        0                                                                     \
    );                                                                        \
                                                                              \
    solidChemistrySolver<Comp, SThermo, GThermo>::                            \
        adddictionaryConstructorToTable<SS<Comp, SThermo, GThermo> >          \
            add##SS##Comp##SThermo##GThermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



#ifdef NoRepository
#   include "solidChemistrySolver.C"
#   include "newSolidChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
