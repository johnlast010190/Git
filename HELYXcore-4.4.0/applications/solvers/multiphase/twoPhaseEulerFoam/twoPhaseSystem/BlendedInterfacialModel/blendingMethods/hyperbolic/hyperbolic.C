/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "hyperbolic.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(hyperbolic, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        hyperbolic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::hyperbolic
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    transitionAlphaScale_
    (
        "transitionAlphaScale",
        dimless,
        dict.lookup("transitionAlphaScale")
    )
{
    forAllConstIter(wordList, phaseNames, iter)
    {
        const word name(IOobject::groupName("maxDispersedAlpha", *iter));

        maxDispersedAlpha_.insert
        (
            *iter,
            dimensionedScalar
            (
                name,
                dimless,
                dict.lookup(name)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::~hyperbolic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase1 - maxDispersedAlpha_[phase1.name()])
            )
        )/2;
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(maxDispersedAlpha_[phase2.name()] - phase2)
            )
        )/2;
}


// ************************************************************************* //
