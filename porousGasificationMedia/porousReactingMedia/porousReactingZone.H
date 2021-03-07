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
    Foam::porousReactingZone

Description
    Porous zone definition based on cell zones.

    Porous zone definition based on cell zones and parameters obtained from a
    control dictionary constructed from the given stream. The orientation of
    the porous region is defined with the same notation as a coordinateSystem,
    but only a Cartesian coordinate system is valid.

    Implemented porosity models:

    Darcy (@e d parameter)
    @f[
        S = - (\mu \, d ) U
    @f]

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The porousReactingZones method porousReactingZones::ddt() mirrors the normal fvm::ddt()
    method, but accounts for the effective volume of the cells.

See Also
    porousReactingZones and coordinateSystems

SourceFiles
    porousReactingZone.C
    porousReactingZoneTemplates.C

    created by Pawel Jan Zuk

\*---------------------------------------------------------------------------*/

#ifndef porousReactingZone_H
#define porousReactingZone_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "IOPtrList.H"
#include "fvMatrix.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class porousReactingZone Declaration
\*---------------------------------------------------------------------------*/

class porousReactingZone
{
    // Private data

        //- Name of this zone
        word name_;

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;

        //- Dictionary containing the parameters
        dictionary dict_;

        //- Cell zone ID
        label cellZoneID_;

        //- Coordinate system used for the zone (Cartesian)
        coordinateSystem coordSys_;

        //- porosity of the zone (0 <= porosity <= 1)
        //  Placeholder for treatment of temporal terms.
        //  Currently unused.
        scalar porosity_;

        //- powerLaw coefficient C0
        scalar C0_;

        //- powerLaw coefficient C1
        scalar C1_;

        //- Darcy coefficient
        dimensionedTensor D_;

        //- Forchheimer coefficient
        dimensionedTensor F_;

        //- porosity field
        volScalarField& porosityF_;

    // Private Member Functions

        //- adjust negative resistance values to be multiplier of max value
        static void adjustNegativeResistance(dimensionedVector& resist);
	
	//const labelList findPorousZone(scalarField& porosityF);


        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            scalarField& Udiag,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            scalarField& Udiag,
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

	//- Df based viscous and intertial resistance
	template<class RhoFieldType>
	void addViscousInertialResistance
	(
	    scalarField& Udiag,
	    vectorField& Usource,
	    const labelList& cells,
	    const scalarField& V,
	    const RhoFieldType& rho,
	    const scalarField& mu,
	    const vectorField& U,
	    tensorField& Df 
	) const;

        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            tensorField& AU,
            const labelList& cells,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            tensorField& AU,
            const labelList& cells,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;


        //- Disallow default bitwise copy construct
        porousReactingZone(
		const porousReactingZone&
		          );

        //- Disallow default bitwise assignment
        void operator=(const porousReactingZone&);

        //- modify time derivative elements
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;

public:

    // Constructors

	//- construct from porosityF 
        porousReactingZone(const fvMesh& mesh,volScalarField& porosityF);

        //- Return clone
        autoPtr<porousReactingZone> clone() const
        {
            notImplemented("autoPtr<porousReactingZone> clone() const");
            return autoPtr<porousReactingZone>(NULL);
        }

    //- Destructor
    virtual ~porousReactingZone()
    {}


    // Member Functions

        // Access

            //- cellZone name
            const word& zoneName() const
            {
                return name_;
            }

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- cellZone number
            label zoneId() const
            {
                return cellZoneID_;
            }

            //- dictionary values used for the porousReactingZone
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return coordinate system
            const coordinateSystem& coordSys() const
            {
                return coordSys_;
            }

            //- Return origin
            const point& origin() const
            {
                return coordSys_.origin();
            }

            //- Return axis
            vector axis() const
            {
                return coordSys_.axis();
            }

            //- Return porosity
            scalar porosity() const
            {
                return porosity_;
            }

            //- Edit access to porosity
            scalar& porosity()
            {
                return porosity_;
            }

	    //- Return porosity indicator field
            volScalarField& porosityF() const
            {
                return porosityF_;
            }

        //- Modify time derivative elements according to porosity
//        template<class Type>
//        void modifyDdt(fvMatrix<Type>&) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
	void addResistance(fvVectorMatrix& UEqn, volTensorField& Df) const;
		
        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal.
        //  Optionally correct the processor BCs of AU.
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        ) const;

        //- Write the porousReactingZone dictionary
        virtual void writeDict(Ostream&, bool subDict = true) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const porousReactingZone&);

	//- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const geometricOneField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        template<class Type>
        tmp<fvMatrix<Type> > ddt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
	);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousReactingZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
