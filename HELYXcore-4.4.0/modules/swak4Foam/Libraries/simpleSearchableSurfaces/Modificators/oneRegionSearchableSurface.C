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

#include "oneRegionSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/Lists/SortableList/SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(oneRegionSearchableSurface, 0);
    addToRunTimeSelectionTable(searchableSurface, oneRegionSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneRegionSearchableSurface::oneRegionSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    wrapperSearchableSurface(io,dict),
    name_(1,word(dict.lookup("regionName")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oneRegionSearchableSurface::~oneRegionSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::wordList& Foam::oneRegionSearchableSurface::regions() const
{
    return name_;
}

#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
Foam::tmp<Foam::pointField>
#else
Foam::pointField
#endif
Foam::oneRegionSearchableSurface::coordinates() const
{
#ifdef FOAM_SEARCHABLE_SURF_USES_TMP
    return tmp<pointField>(new pointField(1,delegate().coordinates()()[0]));
#else
    return pointField(1,delegate().coordinates()[0]);
#endif
}

void Foam::oneRegionSearchableSurface::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region,
    const bool threaded /*= false*/
) const
{
    region.setSize(info.size());
    region=0;
}




// ************************************************************************* //
