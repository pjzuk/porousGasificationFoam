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

\*---------------------------------------------------------------------------*/

#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(heterogeneousAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(heterogeneousAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousAbsorptionEmissionModel::heterogeneousAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousAbsorptionEmissionModel::~heterogeneousAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::a(const label bandI) const
{
    return aDisp(bandI) + aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::aCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::aDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "aDisp",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::as(const label bandI) const
{
    return asCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::asCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "asCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::borderAs(const label bandI) const
{
    return borderAsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::borderAsCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "borderAsCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::e(const label bandI) const
{
    return eDisp(bandI) + eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::eCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::eDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "eDisp",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::E(const label bandI) const
{
    return EDisp(bandI) + ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::es(const label bandI) const
{
    return esCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::esCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "esCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::borderEs(const label bandI) const
{
    return borderEsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::borderEsCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "borderEsCont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::ECont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ECont",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::EDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "EDisp",
                mesh_.time().timeName(),
                mesh_, // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

Foam::dimensionedScalar
Foam::radiation::heterogeneousAbsorptionEmissionModel::borderL(const label bandI) const
{
    return  dimensionedScalar("borderL", dimLength, 0.0);
}

Foam::label Foam::radiation::heterogeneousAbsorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::radiation::heterogeneousAbsorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::radiation::heterogeneousAbsorptionEmissionModel::isGrey() const
{
    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousAbsorptionEmissionModel::addIntensity
(
    const label rayI,
    const volScalarField& ILambda
) const
{
    return ILambda;
}


void Foam::radiation::heterogeneousAbsorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a.internalField() = this->a();
    aj[0].internalField() =  a.internalField();
}


// ************************************************************************* //
