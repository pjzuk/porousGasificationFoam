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

\*---------------------------------------------------------------------------*/

#include "heterogeneousMeanTemp.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

using namespace Foam::radiation;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(heterogeneousMeanTemp, 0);

        addToRunTimeSelectionTable
        (
            heterogeneousRadiationModel,
            heterogeneousMeanTemp,
            porosity
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousMeanTemp::heterogeneousMeanTemp(const volScalarField& T, const volScalarField& porosityF, const volScalarField& surfF, const volScalarField& Ts)
:
    heterogeneousRadiationModel(T, porosityF,  surfF, Ts),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    borderAs_
    (
        IOobject
        (
            "borderAs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("borderAs", dimless/dimLength, 0.0)
    ),
    borderL_
    (
        dimensionedScalar("borderL", dimLength, 0.0)
    ),
    porosityF_
    (
        porosityF
    ),
    surfFI_
    (
        surfF
    ),
    surfF_
    (
        IOobject
        (
            "surfF_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    solidSh_
    (
        IOobject
        (
            "radiationSolidSh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("radiationSolidSh", dimMass/dimLength/pow3(dimTime), 0.0)
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousMeanTemp::~heterogeneousMeanTemp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::heterogeneousMeanTemp::read()
{
    if (heterogeneousRadiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::heterogeneousMeanTemp::calculate()
{
    borderAs_ = absorptionEmission_->borderAs();
    borderL_  = absorptionEmission_->borderL();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    volScalarField surfV = surfF_*0;
    surfV.internalField() = borderL_*pow(mesh_.V(),2./3.)*surfFI_.internalField();
    surfF_.internalField() = borderL_/pow(mesh_.V(),1./3.)*surfFI_.internalField();
    scalar totSur = gSum(surfV);

    scalar totVol = 0;
    forAll(porosityF_,cellI)
    {
        if (porosityF_[cellI] < (1.0 - pow(10.0,-8.0)))  //this is an ad hoc threshold
        {
            totVol += mesh_.V()[cellI];
        }
    }
    reduce(totVol, sumOp<scalar>());

    Info << "Radiation active volume to porous media volume ratio: " << totSur/max(totVol,SMALL) << endl;

    dimensionedScalar boundaryMeanTemp("boundaryMeanTemp",dimless,0.0);
    dimensionedScalar boundarySurface("boundarySurface",dimless,0.0);

    // Calculate radiative heat flux on boundaries.
    bool thereIsWall = false;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "wall")
        {
            boundaryMeanTemp += gSum(T_.boundaryField()[patchI]*mesh_.magSf().boundaryField()[patchI]);
            boundarySurface += gSum(mesh_.magSf().boundaryField()[patchI]);
            thereIsWall = true;
        }
    }

    Info << "Radiation mean wall temperature [K] " << boundaryMeanTemp.value()/boundarySurface.value() << endl;

    if (not thereIsWall)
    {
        FatalErrorIn("Foam::radiation::heterogeneousMeanTemp")
            << "there is no patch of type: wall "
            << nl << "cannot calculate radiation source."
            << exit(FatalError);
    }
    // Solve G transport equation


    // eqZx2uHGn016
    volScalarField solidRadiation = borderAs_*radiation::sigmaSB*pow4(Ts_);

    scalar radiationEnergy = (radiation::sigmaSB*pow4(boundaryMeanTemp/boundarySurface)).value();

    forAll(G_,cellI)
    {
        G_[cellI] = radiationEnergy;
        if (surfF_[cellI] != 0)
        {
            solidSh_[cellI] = 4.0*(G_[cellI]*borderAs_[cellI] - solidRadiation[cellI])*surfF_[cellI];
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::heterogeneousMeanTemp::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Rp",
                radiation::sigmaSB.dimensions()/dimLength,
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::heterogeneousMeanTemp::Ru() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Ru", dimMass/dimLength/pow3(dimTime), 0.0
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousMeanTemp::solidSh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "tSolidSh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidSh_
        )
    );
}


// ************************************************************************* //
