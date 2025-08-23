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
    (c) 2014-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "linear.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianBlendingMethods
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        eulerianBlendingMethod,
        linear,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianBlendingMethods::linear::linear
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    eulerianBlendingMethod(dict)
{
    forAllConstIter(wordList, phaseNames, iter)
    {
        const word memberFull("minFullyContinuousAlpha");
        const word nameFull(IOobject::groupName(memberFull, *iter));

        minFullyContinuousAlpha_.insert
        (
            *iter,
            dimensionedScalar
            (
                nameFull,
                dimless,
                (!dict.found(nameFull) && dict.found(memberFull))
              ? dict.lookup(memberFull)
              : dict.lookup(nameFull)
            )
        );

        const word memberPart("minPartlyContinuousAlpha");
        const word namePart
        (
            IOobject::groupName(memberPart, *iter)
        );

        minPartlyContinuousAlpha_.insert
        (
            *iter,
            dimensionedScalar
            (
                namePart,
                dimless,
                (!dict.found(namePart) && dict.found(memberPart))
              ? dict.lookup(memberPart)
              : dict.lookup(namePart)
            )
        );

        if
        (
            minFullyContinuousAlpha_[*iter]
          < minPartlyContinuousAlpha_[*iter]
        )
        {
            FatalErrorInFunction
                << "The supplied fully continuous volume fraction for "
                << *iter
                << " is less than the partly continuous value."
                << endl << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianBlendingMethods::linear::~linear()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianBlendingMethods::linear::f1
(
    const eulerianPhaseModel& phase1,
    const eulerianPhaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase2.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase2.name()]);

    return
        min
        (
            max
            (
                (phase2.volFrac() - minPartAlpha)
               /(minFullAlpha - minPartAlpha + SMALL),
                scalar(0)
            ),
            scalar(1)
        );
}


Foam::tmp<Foam::volScalarField> Foam::eulerianBlendingMethods::linear::f2
(
    const eulerianPhaseModel& phase1,
    const eulerianPhaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase1.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase1.name()]);

    return
        min
        (
            max
            (
                (phase1.volFrac() - minPartAlpha)
               /(minFullAlpha - minPartAlpha + SMALL),
                scalar(0)
            ),
            scalar(1)
        );
}


// ************************************************************************* //
