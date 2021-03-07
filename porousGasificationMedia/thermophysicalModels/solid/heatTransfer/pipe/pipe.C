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

#include "pipe.H"
#include "foamTime.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
      

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pipeCONV, 0);
addToRunTimeSelectionTable(heatTransferModel, pipeCONV, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pipeCONV::pipeCONV
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    heatTransferModel(por,por0),
    hCoeff_(1.0),
    pipeRadius_(1.0),
    Up_(db().lookupObject<volVectorField>("U")),
    rhop_(db().lookupObject<volScalarField>("rho")),
    alphap_(db().lookupObject<volScalarField>("alpha")),
    mup_(db().lookupObject<volScalarField>("mu")),
    thermop_(db().lookupObject<basicThermo>("thermophysicalProperties"))
{
   read();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<pipeCONV> pipeCONV::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<pipeCONV>
    (
        new pipeCONV( por,por0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> pipeCONV::CONV() const
{
// eqZx2uHGn007
    Foam::tmp<Foam::volScalarField> CONVloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONVloc",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", dimEnergy/dimTime/dimTemperature/dimVolume, 0.0
            )
        )
    );
    
    volScalarField& Cp = thermop_.Cp()();
    forAll (CONVloc_(),cellI)
    {
	    CONVloc_()[cellI] = pow(por()[cellI],0.5)*pow(por0()[cellI],0.5)*2.0/pipeRadius_*
              (3.66)
                *Cp[cellI]*alphap_[cellI]*rhop_[cellI]/pipeRadius_;  //eqZx2uHGn019 eqZx2uHGn020 
    }

    return CONVloc_;
}

bool pipeCONV::read()
{

	IOdictionary dict
        (
            IOobject
            (
                "heatTransferProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    const dictionary& params = dict.subDict("Parameters");

    params.lookup("h") >> hCoeff_;
    params.lookup("pipeRadius") >> pipeRadius_;

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
