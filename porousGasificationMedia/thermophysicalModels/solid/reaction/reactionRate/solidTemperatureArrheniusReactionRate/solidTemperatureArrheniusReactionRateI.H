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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidTemperatureArrheniusReactionRate::solidTemperatureArrheniusReactionRate
(
    const scalar A,
    const scalar Ta,
    const scalar Tcrit
)
:
    A_(A),
    Ta_(Ta),
    Tcrit_(Tcrit)
{}


inline Foam::solidTemperatureArrheniusReactionRate::solidTemperatureArrheniusReactionRate
(
    const speciesTable&,
    Istream& is
)
{
    is.readBegin("solidTemperatureArrheniusReaction(Istream&)");
    A_ = scalar(readScalar(is));
    Ta_ = scalar(readScalar(is));
    Tcrit_ = scalar(readScalar(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidTemperatureArrheniusReactionRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField&
) const

// eqId (5)
{
    scalar ak = A_;

    if (T < Tcrit_)
    {
        ak *= 0.0;
    }
    else
    {
        ak *= T*exp(-Ta_/T);
    }

  return ak;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidTemperatureArrheniusReactionRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.Ta_ << token::SPACE << arr.Tcrit_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
