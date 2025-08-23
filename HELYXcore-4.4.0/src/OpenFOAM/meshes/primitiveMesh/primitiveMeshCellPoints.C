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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellPoints() const
{
    HELYX_ASSERT(!cpPtr_)
    {
        FatalErrorInFunction << "Should not re-calculate cell points." << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "primitiveMesh::cellPoints() : "
            << "calculating cellPoints" << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // Invert pointCells
    cpPtr_ = new labelListList(nCells());
    invertManyToMany(nCells(), pointCells(), *cpPtr_);
}


const Foam::labelList& Foam::primitiveMesh::cellPoints
(
    const label celli,
    DynamicList<label>& storage
) const
{
    if (hasCellPoints())
    {
        return cellPoints()[celli];
    }
    else
    {
        const faceList& fcs = faces();
        const labelList& cFaces = cells()[celli];

        labelSet_.clear();

        forAll(cFaces, i)
        {
            const labelList& f = fcs[cFaces[i]];

            forAll(f, fp)
            {
                labelSet_.insert(f[fp]);
            }
        }

        storage.clear();
        if (labelSet_.size() > storage.capacity())
        {
            storage.setCapacity(labelSet_.size());
        }

        forAllConstIter(labelHashSet, labelSet_, iter)
        {
            storage.append(iter.key());
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::cellPoints(const label celli) const
{
    return cellPoints(celli, labels_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
