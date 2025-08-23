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
    2009, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "translateSearchableSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(translateSearchableSurface, 0);
addToRunTimeSelectionTable(searchableSurface, translateSearchableSurface, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translateSearchableSurface::translateSearchableSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    transformationSearchableSurface(io,dict),
    translation_(dict.lookup("translation"))
{
#ifdef FOAM_SEARCHABLE_SURF_HAS_BOUND_METHOD
    bounds()=boundBox(
        transform(delegate().bounds().min()),
        transform(delegate().bounds().max())
    );
#endif
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::translateSearchableSurface::~translateSearchableSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::translateSearchableSurface::transform(const point &p) const
{
    return p+translation_;
}

Foam::point Foam::translateSearchableSurface::inverseTransform(const point &p) const
{
    return p-translation_;
}

#ifdef FOAM_SEARCHABLE_SURF_NEEDS_BOUNDING_SPHERES
void Foam::translateSearchableSurface::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    delegate().boundingSpheres(
        centres,
        radiusSqr
    );
    forAll(centres,i) {
        centres[i]=centres[i]+translation_;
    }
}
#endif

// ************************************************************************* //
