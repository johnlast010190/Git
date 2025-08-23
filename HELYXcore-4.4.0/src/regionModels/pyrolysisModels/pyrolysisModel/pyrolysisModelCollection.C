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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "pyrolysisModels/pyrolysisModel/pyrolysisModelCollection.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(pyrolysisModelCollection, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisModelCollection::pyrolysisModelCollection
(
    const fvMesh& mesh
)
:
    PtrList<pyrolysisModel>()
{
    IOdictionary pyrolysisZonesDict
    (
        IOobject
        (
            "pyrolysisZones",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const wordList regions(pyrolysisZonesDict.toc());

    setSize(regions.size());
    forAll(regions, regioni)
    {
        const dictionary& regionDict =
            pyrolysisZonesDict.subDict(regions[regioni]);
        set(regioni, pyrolysisModel::New(mesh, regionDict, regions[regioni]));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisModelCollection::~pyrolysisModelCollection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionModels::pyrolysisModelCollection::evolveRegion()
{
    forAll(*this, i)
    {
        this->operator[](i).evolveRegion();
    }
}


void Foam::regionModels::pyrolysisModelCollection::evolve()
{
    forAll(*this, i)
    {
        pyrolysisModel& pyrolysis = this->operator[](i);

        if (pyrolysis.active())
        {
            if (pyrolysis.primaryMesh().changing())
            {
                FatalErrorInFunction
                    << "Currently not possible to apply "
                    << pyrolysis.modelName()
                    << " model to moving mesh cases" << nl<< abort(FatalError);
            }

            // Increment the region equations up to the new time level
            pyrolysis.evolveRegion();

            // Provide some feedback
            if (pyrolysis.infoOutput())
            {
                Info<< incrIndent;
                pyrolysis.info();
                Info<< endl << decrIndent;
            }
        }
    }
}


void Foam::regionModels::pyrolysisModelCollection::info()
{
    forAll(*this, i)
    {
        this->operator[](i).info();
    }
}


Foam::scalar Foam::regionModels::pyrolysisModelCollection::maxDiff() const
{
    scalar maxDiff = 0.0;
    forAll(*this, i)
    {
        maxDiff = max(maxDiff, this->operator[](i).maxDiff());
    }

    return maxDiff;
}


Foam::scalar
Foam::regionModels::pyrolysisModelCollection::solidRegionDiffNo() const
{
    scalar totalDiNum = GREAT;
    forAll(*this, i)
    {
        totalDiNum = min(totalDiNum, this->operator[](i).solidRegionDiffNo());
    }
    if (totalDiNum == -GREAT)
    {
        totalDiNum = SMALL;
    }

    return totalDiNum;
}


// ************************************************************************* //
