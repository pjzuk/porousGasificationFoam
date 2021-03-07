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

// adjust negative resistance values to be multiplier of max value
void Foam::porousReactingZone::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorIn
        (
            "Foam::porousReactingZone::porousReactingZone::adjustNegativeResistance"
            "(dimensionedVector&)"
        )   << "negative resistances! " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt=0; cmpt < vector::nComponents; ++cmpt)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


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
    C0_(0.),
    C1_(0.),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero),
    porosityF_(porosityF)
{
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



void Foam::porousReactingZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                AU,
                cells,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
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

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        os << indent << "powerLaw";
        dictPtr->write(os);
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
