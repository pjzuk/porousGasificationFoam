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

\*----------------------------------------------------------------------------*/

#include "porousReactingZone.H"
#include "List.H"
#include "label.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "foamTime.H"
#include "volFields.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::porousReactingZone::porousReactingZone
(
    const fvMesh& mesh,
    volScalarField& porosityF
)
:
    name_(),
    mesh_(mesh),
    dict_(),
    cellZoneID_(),
    coordSys_(),
    porosity_(0.),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    f_(0.),
    porosityF_(porosityF)
{

    IOobject porosityPropertiesHeader
    (
        "porosityProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ
    );

    word modelType;
    if (porosityPropertiesHeader.headerOk())
    {
        IOdictionary dict
        (
            IOobject
            (
                "porosityProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        f_ = dict.lookupOrDefault("forchheimerCoeff",0.);
        Info << "Forchheimer coefficient f = " << f_  << " specified" << nl 
             << "Darcy resistance term will be modified by + f*rho*mag(u)*sqrt(3.)*Df/|Df|" 
             << nl << endl;

    }
    else
    {
        Info << "no Forchheimer coefficient specified" << nl 
             << "Darcy resistance term only. For the Forchheimer term create porosityProperties dictionary." 
             << nl << endl;
    }

}

//*****************************************************************************************

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousReactingZone::addResistance
(
	fvVectorMatrix& UEqn,
	volTensorField& Df
)
const
{
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    
	label n=0;
	forAll(porosityF_,celli)
		if(porosityF_[celli] >= 0 && porosityF_[celli] < 1)
			n++;
	labelList cells(n);
	label i=0;
	forAll(porosityF_,celli)
	{
		if(porosityF_[celli] >= 0 && porosityF_[celli] < 1)
		{
			cells[i]=celli;
			i++;
		}
	}
    
    addViscousInertialResistance
	        (
	            Udiag,
	            Usource,
	            cells,
	            V,
	            mesh_.lookupObject<volScalarField>("rho"),
	            mesh_.lookupObject<volScalarField>("mu"),
	            U,
	            Df
	        );
}


void Foam::porousReactingZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name") << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity") << porosity() << token::END_STATEMENT << nl;
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        os << indent << "Darcy";
        dictPtr->write(os);
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porousReactingZone& pZone)
{
    pZone.writeDict(os);
    return os;
}

// ************************************************************************* //
