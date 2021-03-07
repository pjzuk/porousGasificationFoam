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

#include "basicHGSSolidThermo.H"
#include "fvMesh.H"
#include "HashTable.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicHGSSolidThermo, 0);
    defineRunTimeSelectionTable(basicHGSSolidThermo, mesh);
    defineRunTimeSelectionTable(basicHGSSolidThermo, dictionary);
    defineRunTimeSelectionTable(basicHGSSolidThermo, gasPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicHGSSolidThermo::basicHGSSolidThermo(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappas",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )

{
}


Foam::basicHGSSolidThermo::basicHGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{
}


Foam::basicHGSSolidThermo::basicHGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const PtrList<volScalarField>& gasPhaseGases
)
:
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicHGSSolidThermo::~basicHGSSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicHGSSolidThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicHGSSolidThermo::T() const
{
    return T_;
}


const Foam::volScalarField& Foam::basicHGSSolidThermo::rho() const
{
    return rho_;
}


Foam::volScalarField& Foam::basicHGSSolidThermo::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::basicHGSSolidThermo::kappa() const
{
    return kappa_;
}


const Foam::volScalarField& Foam::basicHGSSolidThermo::sigmaS() const
{
    return sigmaS_;
}


const Foam::volScalarField& Foam::basicHGSSolidThermo::emissivity() const
{
    return emissivity_;
}


const Foam::volScalarField&  Foam::basicHGSSolidThermo::K() const
{
    notImplemented("basicHGSSolidThermo::K()");
    return volScalarField::null();
}


const Foam::volSymmTensorField& Foam::basicHGSSolidThermo::directionalK() const
{
    notImplemented("basicHGSSolidThermo::directionalK()");
    return const_cast<volSymmTensorField&>(volSymmTensorField::null());
}


Foam::basicSolidMixture& Foam::basicHGSSolidThermo::composition()
{
    notImplemented("basicHGSSolidThermo::composition()");
    return *reinterpret_cast<basicSolidMixture*>(0);
}


const Foam::basicSolidMixture& Foam::basicHGSSolidThermo::composition() const
{
    notImplemented("basicHGSSolidThermo::composition() const");
    return *reinterpret_cast<const basicSolidMixture*>(0);
}


Foam::tmp<Foam::volScalarField> Foam::basicHGSSolidThermo::hs() const
{
    notImplemented("basicHGSSolidThermo::hs()");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicHGSSolidThermo::hs(const label patchI)
const
{
    notImplemented("basicHGSSolidThermo::hs(const label)");
    return scalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicHGSSolidThermo::K
(
    const label patchI
)const
{
    notImplemented("basicHGSSolidThermo::K(const label)");
    return scalarField::null();
}


Foam::tmp<Foam::symmTensorField> Foam::basicHGSSolidThermo::directionalK
(
    const label
)const
{
    notImplemented("basicHGSSolidThermo::directionalK(const label)");
    return symmTensorField::null();
}


bool Foam::basicHGSSolidThermo::read()
{
    return regIOobject::read();
}


bool Foam::basicHGSSolidThermo::writeData(Ostream& os) const
{
    return true;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const basicHGSSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
