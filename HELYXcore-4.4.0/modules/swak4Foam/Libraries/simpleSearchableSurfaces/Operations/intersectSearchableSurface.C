/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "intersectSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(intersectSearchableSurface, 0);
addToRunTimeSelectionTable(searchableSurface, intersectSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::intersectSearchableSurface::intersectSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    binaryOperationSearchableSurface(io,dict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::intersectSearchableSurface::~intersectSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::intersectSearchableSurface::decidePoint(
        const hitWhom who,
        const bool inA,
        const bool inB
) const
{
    if (
        who==BOTH
    ) {
        return true;
    } else if (
        (inA && who==HITSB)
        ||
        (inB && who==HITSA)
    ) {
        return true;
    }

    return false;
}

void Foam::intersectSearchableSurface::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType,
    const bool threaded /*= false*/
) const
{
    List<volumeType> inA;
    List<volumeType> inB;

    a().getVolumeType(points,inA);
    b().getVolumeType(points,inB);

    volType.setSize(points.size());

    forAll(volType,i) {
        if (inA[i]==INSIDE && inB[i]==INSIDE) {
            volType[i]=INSIDE;
        } else if (inA[i]==UNKNOWN || inB[i]==UNKNOWN) {
            volType[i]=UNKNOWN;
        } else {
            volType[i]=OUTSIDE;
        }
    }
}

// ************************************************************************* //
