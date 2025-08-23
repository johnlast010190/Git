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

#include "invertSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(invertSearchableSurface, 0);
    addToRunTimeSelectionTable(searchableSurface, invertSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::invertSearchableSurface::invertSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    wrapperSearchableSurface(io,dict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::invertSearchableSurface::~invertSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::invertSearchableSurface::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType,
    const bool threaded /*= false*/
) const
{
    delegate().getVolumeType
        (
            points,
            volType
        );

    forAll(volType,i) {
        if (volType[i]==INSIDE) {
            volType[i]=OUTSIDE;
        } else if (volType[i]==OUTSIDE) {
            volType[i]=INSIDE;
        }
    }
}

void Foam::invertSearchableSurface::getNormal
(
    const List<pointIndexHit>& hits,
    vectorField& normal
) const
{
    delegate().getNormal(
        hits,
        normal
    );

    forAll(normal,i) {
        normal[i]=-normal[i];
    }
}

// ************************************************************************* //
