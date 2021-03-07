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

Class
    Foam::solidChemistryModel

Description
    Chemistry model for solid thermodynamics

SourceFiles
    solidChemistryModelI.H
    solidChemistryModel.C
    newChemistrySolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidChemistryModel_H
#define solidChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "basicHGSSolidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class solidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class solidChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        solidChemistryModel(const solidChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidChemistryModel&);


protected:

    // Protected data

        //- Solid thermo package
        autoPtr<basicHGSSolidThermo> solidThermo_;


public:

    //- Runtime type information
    TypeName("solidChemistryModel");


    //- Declare run-time constructor selection tables
/*    declareRunTimeSelectionTable
    (
        autoPtr,
        solidChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
	    const objectRegistry& obj,
            const word& compTypeName,
            const word& solidThermoTypeName
        ),
        (mesh, obj, compTypeName, solidThermoTypeName)
    );
*/
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidChemistryModel,
        fvMeshGasPhaseGases,
        (
            const fvMesh& mesh,
	    const objectRegistry& obj,
            const word& compTypeName,
            const word& solidThermoTypeName,
	    PtrList<volScalarField>& gasPhaseGases
        ),
        (mesh, obj, compTypeName, solidThermoTypeName, gasPhaseGases)
    );

    // Constructors

        //- Construct from mesh and thermo type name
        solidChemistryModel
        (
            const fvMesh& mesh,
	    const objectRegistry& obj,
	    const word& compTypeName,
            const word& solidThermoTypeName
        );

	solidChemistryModel
	(
	    const fvMesh& mesh,
	    const objectRegistry& obj,
	    const word& compTypeName,
	    const word& solidThermoTypeName,
	    PtrList<volScalarField>& gasPhaseGases
	);

    //- Selectors
/*
    static autoPtr<solidChemistryModel> New(const fvMesh&, const objectRegistry&);

    static autoPtr<solidChemistryModel> New(const fvMesh& mesh)
    {
        return New(mesh, mesh);
    }
*/
    static autoPtr<solidChemistryModel> New(const fvMesh& mesh, const objectRegistry&, PtrList<volScalarField>& gasPhaseGases);

    static autoPtr<solidChemistryModel> New(const fvMesh& mesh, PtrList<volScalarField>& gasPhaseGases)
    {
        return New(mesh, mesh, gasPhaseGases);
    }


    //- Destructor
    virtual ~solidChemistryModel();


    // Member Functions

        //- Return access to the solid thermo package
        inline basicHGSSolidThermo& solidThermo();

        //- Return const access to the solid thermo package
        inline const basicHGSSolidThermo& solidThermo() const;

        //- Return total gases mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRg() const = 0;

        //- Return total solids mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRs() const = 0;

        //- Return chemical source terms for solids [kg/m3/s]
        virtual tmp<volScalarField> RRs(const label i) const = 0;

        //- Return chemical source terms for gases [kg/m3/s]
        virtual tmp<volScalarField> RRg(const label i) const = 0;

        //- Return sensible enthalpy for gas i [J/Kg]
        virtual tmp<volScalarField> gasHs
        (
            const volScalarField& T,
            const label i
        ) const = 0;

        //- Return specie Table for gases
        virtual const speciesTable& gasTable() const = 0;

        //- Set reacting status of cell, cellI
        virtual void setCellReacting(const label cellI, const bool active) = 0;

        //- Return source for porosity equation [1/s]
        virtual tmp<volScalarField> RRpor(const volScalarField T) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
