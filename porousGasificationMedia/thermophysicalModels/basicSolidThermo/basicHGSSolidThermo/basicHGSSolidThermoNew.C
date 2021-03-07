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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicHGSSolidThermo> Foam::basicHGSSolidThermo::New
(
    const fvMesh& mesh
)
{
    if (debug)
    {
        Info<< "basicHGSSolidThermo::New(const fvMesh&): "
            << "constructing basicHGSSolidThermo"
            << endl;
    }

    const word thermoType
    (
        IOdictionary
        (
            IOobject
            (
                "solidThermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("thermoType")
    );

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(thermoType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicHGSSolidThermo::New(const fvMesh&)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << meshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicHGSSolidThermo>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::basicHGSSolidThermo> Foam::basicHGSSolidThermo::New
(
    const fvMesh& mesh, const dictionary& dict
)
{
    if (debug)
    {
        Info<< "basicHGSSolidThermo::New(const fvMesh&, const dictionary&): "
            << "constructing basicHGSSolidThermo"
            << endl;
    }

    const word thermoType = dict.lookup("thermoType");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermoType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicHGSSolidThermo::New(const fvMesh&, const dictionary&)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicHGSSolidThermo>(cstrIter()(mesh, dict));
}


Foam::autoPtr<Foam::basicHGSSolidThermo> Foam::basicHGSSolidThermo::New
(
    const fvMesh& mesh, const dictionary& dict, const PtrList<volScalarField>& gasPhaseGases
)
{
    if (debug)
    {
        Info<< "basicHGSSolidThermo::New(const fvMesh&, const dictionary&, PtrList<volScalarField>& gasPhaseGases): "
            << "constructing basicHGSSolidThermo"
            << endl;
    }

    const word thermoType = dict.lookup("thermoType");

    gasPhaseConstructorTable::iterator cstrIter =
        gasPhaseConstructorTablePtr_->find(thermoType);

    if (cstrIter == gasPhaseConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicHGSSolidThermo::New(const fvMesh&, const dictionary&, PtrList<volScalarField>& gasPhaseGases)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << gasPhaseConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicHGSSolidThermo>(cstrIter()(mesh, dict, gasPhaseGases));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
