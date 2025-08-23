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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2009, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "searchableBoxWithDirections.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableBoxWithDirections, 0);
addToRunTimeSelectionTable(searchableSurface, searchableBoxWithDirections, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableBoxWithDirections::searchableBoxWithDirections
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableBox(
        io,
        dict
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableBoxWithDirections::~searchableBoxWithDirections()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableBoxWithDirections::regions() const
{
    if (regions_.size() == 0)
    {
        regions_.setSize(6);
        regions_[0] = "xMin";
        regions_[1] = "xMax";
        regions_[2] = "yMin";
        regions_[3] = "yMax";
        regions_[4] = "zMin";
        regions_[5] = "zMax";
    }
    return regions_;
}


void Foam::searchableBoxWithDirections::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region,
    const bool threaded /*= false*/
) const
{
    region.setSize(info.size());
    forAll(info,i) {
        region[i] = info[i].index();
    }
}

// ************************************************************************* //
