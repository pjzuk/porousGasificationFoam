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

#include "solidChemistrySolver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class SThermoType, class GThermoType>
Foam::autoPtr<Foam::solidChemistrySolver<CompType, SThermoType, GThermoType> >
Foam::solidChemistrySolver<CompType, SThermoType, GThermoType>::New
(
    ODESolidHeterogeneousChemistryModel<CompType, SThermoType, GThermoType>& model,
    const word& compTypeName,
    const word& thermoTypeName
)
{
    word modelName(model.lookup("solidChemistrySolver"));

    word chemistrySolverType =
        modelName + '<' + "solidChemistryModel" + ',' + thermoTypeName + '>';

    Info<< "Selecting solidChemistrySolver " << modelName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistrySolverType);

    Info << modelName << endl;
    Info << compTypeName << endl;
    Info << thermoTypeName << endl;
    Info << chemistrySolverType << endl;

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        wordList models = dictionaryConstructorTablePtr_->sortedToc();

	Info << models << endl;

        forAll(models, i)
        {
            models[i] = models[i].replace
            (
                '<' + compTypeName + ',' + thermoTypeName + '>',
                ""
            );
        }

        FatalErrorIn
        (
            "solidChemistrySolver::New"
            "("
                "const ODEChemistryModel&, "
                "const word&, "
                "const word&"
            ")"
        )   << "Unknown solidChemistrySolver type " << modelName
            << nl << nl << "Valid solidChemistrySolver types are:" << nl
            << models << nl << exit(FatalError);
    }

    return autoPtr<solidChemistrySolver<CompType, SThermoType, GThermoType> >
        (cstrIter()(model, modelName));
}


// ************************************************************************* //
