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
    (c) 1991-2010 OpenCFD Ltd.

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "swakRegistryProxySet.H"
#include "db/dictionary/dictionary.H"
#include "fields/volFields/volFields.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvMesh.H"

#include "repositories/SetsRepository.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swakRegistryProxySet, 0);
    addToRunTimeSelectionTable(sampledSet, swakRegistryProxySet, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void Foam::swakRegistryProxySet::createGeometry()
// {
//     if (debug)
//     {
//         Pout<< "swakRegistryProxySet::createGeometry() - doing nothing"
//             << endl;
//     }
// }

Foam::sampledSet &Foam::swakRegistryProxySet::realSet()
{
    return SetsRepository::getRepository(
        mesh()
    ).getSet(
        setName_,
        static_cast<const fvMesh&>(mesh())
    );
}

const Foam::sampledSet &Foam::swakRegistryProxySet::realSet() const
{
    return SetsRepository::getRepository(
        mesh()
    ).getSet(
        setName_,
        static_cast<const fvMesh&>(mesh())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swakRegistryProxySet::swakRegistryProxySet
(
    const word& name,
    const polyMesh& mesh,
#ifdef FOAM_MESHSEARCH_CONST_SAMPLEDSET
    const meshSearch& search,
#else
    meshSearch& search,
#endif
    const dictionary& dict
)
:
    sampledSet(
        name,
        mesh,
        search,
        dict
    ),
    setName_(dict.lookup("setName"))
{
    setSamples(
        realSet(),
        realSet().cells(),
        realSet().faces(),
        realSet().segments(),
        realSet().curveDist()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swakRegistryProxySet::~swakRegistryProxySet()
{}

#ifdef FOAM_SAMPLEDSET_NEEDS_REFPOINT
Foam::point Foam::swakRegistryProxySet::getRefPoint (const List< point > &pl) const
{
    return realSet().getRefPoint(pl);
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
