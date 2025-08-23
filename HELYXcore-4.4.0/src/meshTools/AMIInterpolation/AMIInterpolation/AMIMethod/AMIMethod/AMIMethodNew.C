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
    (c) 2013-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::AMIMethod> Foam::AMIMethod::New
(
    const word& methodName,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const scalar& cosMatchAngle,
    const bool reverseTarget,
    const bool requireMatch
)
{
    if (debug)
    {
        Info<< "Selecting AMIMethod " << methodName << endl;
    }

    const auto ctor =
        ctorTableLookup
        (
            "AMIMethod type",
            componentsConstructorTable_(),
            methodName
        );
    return autoPtr<AMIMethod>
    (
        ctor
        (
            srcPatch,
            tgtPatch,
            srcMagSf,
            tgtMagSf,
            triMode,
            cosMatchAngle,
            reverseTarget,
            requireMatch
        )
    );
}


// ************************************************************************* //
