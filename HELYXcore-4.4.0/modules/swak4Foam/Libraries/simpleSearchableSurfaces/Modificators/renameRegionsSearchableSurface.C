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

#include "renameRegionsSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "primitives/Pair/Pair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(renameRegionsSearchableSurface, 0);
    addToRunTimeSelectionTable(searchableSurface, renameRegionsSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::renameRegionsSearchableSurface::renameRegionsSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    wrapperSearchableSurface(io,dict)
{
    regions_=delegate().regions();
    List<Pair<word>> replace(dict.lookup("replacementPairs"));

    if (debug) {
        Info<< "renameRegionsSearchableSurface::renameRegionsSearchableSurface"
            << endl;
        Info<< "Old regions: " << regions_ << endl;
        Info<< "Replacements: " << replace << endl;
    }
    forAll(replace,i) {
        if (replace[i].second()=="_") {
            replace[i].second()="";
        }
        forAll(regions_,j) {
            regions_[j].replaceAll
                (
                    replace[i].first(),
                    replace[i].second()
                );
        }
    }
    if (debug) {
        Info<< "New regions: " << regions_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::renameRegionsSearchableSurface::~renameRegionsSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::wordList& Foam::renameRegionsSearchableSurface::regions() const
{
    return regions_;
}



// ************************************************************************* //
