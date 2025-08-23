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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyTopoChangeMap/cellMapper/cellMapper.H"
#include "include/demandDrivenData.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellMapper::calcAddressing() const
{
    if
    (
        directAddrPtr_
     || interpolationAddrPtr_
     || weightsPtr_
     || insertedCellLabelsPtr_
    )
    {
        FatalErrorInFunction
            << "Addressing already calculated."
            << abort(FatalError);
    }

    if (direct())
    {
        // Direct addressing, no weights

        directAddrPtr_ = new labelList(map_.cellMap());
        labelList& directAddr = *directAddrPtr_;

        // Not necessary to resize the list as there are no retired cells
        // directAddr.setSize(mesh_.nCells());

        insertedCellLabelsPtr_ = new labelList(mesh_.nCells());
        labelList& insertedCells = *insertedCellLabelsPtr_;

        label nInsertedCells = 0;

        forAll(directAddr, celli)
        {
            if (directAddr[celli] < 0)
            {
                // Found inserted cell
                directAddr[celli] = 0;
                insertedCells[nInsertedCells] = celli;
                nInsertedCells++;
            }
        }

        insertedCells.setSize(nInsertedCells);
    }
    else
    {
        // Interpolative addressing

        interpolationAddrPtr_ = new labelListList(mesh_.nCells());
        labelListList& addr = *interpolationAddrPtr_;

        weightsPtr_ = new scalarListList(mesh_.nCells());
        scalarListList& w = *weightsPtr_;

        const List<objectMap>& cfp = map_.cellsFromPointsMap();

        forAll(cfp, cfpI)
        {
            // Get addressing
            const labelList& mo = cfp[cfpI].masterObjects();

            label celli = cfp[cfpI].index();

            if (addr[celli].size())
            {
                FatalErrorInFunction
                    << "Master cell " << celli
                    << " mapped from point cells " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[celli] = mo;
            w[celli] = scalarList(mo.size(), 1.0/mo.size());
        }

        const List<objectMap>& cfe = map_.cellsFromEdgesMap();

        forAll(cfe, cfeI)
        {
            // Get addressing
            const labelList& mo = cfe[cfeI].masterObjects();

            label celli = cfe[cfeI].index();

            if (addr[celli].size())
            {
                FatalErrorInFunction
                    << "Master cell " << celli
                    << " mapped from edge cells " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[celli] = mo;
            w[celli] = scalarList(mo.size(), 1.0/mo.size());
        }

        const List<objectMap>& cff = map_.cellsFromFacesMap();

        forAll(cff, cffI)
        {
            // Get addressing
            const labelList& mo = cff[cffI].masterObjects();

            label celli = cff[cffI].index();

            if (addr[celli].size())
            {
                FatalErrorInFunction
                    << "Master cell " << celli
                    << " mapped from face cells " << mo
                    << " already destination of mapping." << abort(FatalError);
            }

            // Map from masters, uniform weights
            addr[celli] = mo;
            w[celli] = scalarList(mo.size(), 1.0/mo.size());
        }

        // Volume conservative mapping if possible

        const List<objectMap>& cfc = map_.cellsFromCellsMap();

        forAll(cfc, cfcI)
        {
            // Get addressing
            const labelList& mo = cfc[cfcI].masterObjects();

            label celli = cfc[cfcI].index();

            if (addr[celli].size())
            {
                FatalErrorInFunction
                    << "Master cell " << celli
                    << " mapped from cell cells " << mo
                    << " already destination of mapping."
                    << abort(FatalError);
            }

            // Map from masters
            addr[celli] = mo;
        }

        if (map_.hasOldCellVolumes())
        {
            // Volume weighted

            const scalarField& V = map_.oldCellVolumes();

            if (V.size() != sizeBeforeMapping())
            {
                FatalErrorInFunction
                    << "cellVolumes size " << V.size()
                    << " is not the old number of cells " << sizeBeforeMapping()
                    << ". Are your cellVolumes already mapped?"
                    << " (new number of cells " << map_.cellMap().size() << ")"
                    << abort(FatalError);
            }

            forAll(cfc, cfcI)
            {
                const labelList& mo = cfc[cfcI].masterObjects();

                label celli = cfc[cfcI].index();

                w[celli].setSize(mo.size());

                if (mo.size())
                {
                    scalar sumV = 0;
                    forAll(mo, ci)
                    {
                        w[celli][ci] = V[mo[ci]];
                        sumV += V[mo[ci]];
                    }
                    if (sumV > VSMALL)
                    {
                        forAll(mo, ci)
                        {
                            w[celli][ci] /= sumV;
                        }
                    }
                    else
                    {
                        // Exception: zero volume. Use uniform mapping
                        w[celli] = scalarList(mo.size(), 1.0/mo.size());
                    }
                }
            }
        }
        else
        {
            // Uniform weighted

            forAll(cfc, cfcI)
            {
                const labelList& mo = cfc[cfcI].masterObjects();

                label celli = cfc[cfcI].index();

                w[celli] = scalarList(mo.size(), 1.0/mo.size());
            }
        }


        // Do mapped faces. Note that can already be set from cellsFromCells
        // so check if addressing size still zero.

        const labelList& cm = map_.cellMap();

        forAll(cm, celli)
        {
            if (cm[celli] > -1 && addr[celli].empty())
            {
                // Mapped from a single cell
                addr[celli] = labelList(1, cm[celli]);
                w[celli] = scalarList(1, 1.0);
            }
        }

        // Grab inserted points (for them the size of addressing is still zero)

        insertedCellLabelsPtr_ = new labelList(mesh_.nCells());
        labelList& insertedCells = *insertedCellLabelsPtr_;

        label nInsertedCells = 0;

        forAll(addr, celli)
        {
            if (addr[celli].empty())
            {
                // Mapped from a dummy cell
                addr[celli] = labelList(1, label(0));
                w[celli] = scalarList(1, 1.0);

                insertedCells[nInsertedCells] = celli;
                nInsertedCells++;
            }
        }

        insertedCells.setSize(nInsertedCells);
    }
}


void Foam::cellMapper::clearOut()
{
    deleteDemandDrivenData(directAddrPtr_);
    deleteDemandDrivenData(interpolationAddrPtr_);
    deleteDemandDrivenData(weightsPtr_);
    deleteDemandDrivenData(insertedCellLabelsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellMapper::cellMapper(const polyTopoChangeMap& map)
:
    mesh_(map.mesh()),
    map_(map),
    insertedCells_(true),
    direct_(false),
    directAddrPtr_(nullptr),
    interpolationAddrPtr_(nullptr),
    weightsPtr_(nullptr),
    insertedCellLabelsPtr_(nullptr)
{
    // Check for possibility of direct mapping
    if
    (
        map_.cellsFromPointsMap().empty()
     && map_.cellsFromEdgesMap().empty()
     && map_.cellsFromFacesMap().empty()
     && map_.cellsFromCellsMap().empty()
    )
    {
        direct_ = true;
    }
    else
    {
        direct_ = false;
    }

    // Check for inserted cells
    if (direct_ && (map_.cellMap().empty() || min(map_.cellMap()) > -1))
    {
        insertedCells_ = false;
    }
    else
    {
        // Need to check all 3 lists to see if there are inserted cells
        // with no owner

        // Make a copy of the cell map, add the entried for cells from points,
        // cells from edges and cells from faces and check for left-overs
        labelList cm(mesh_.nCells(), -1);

        const List<objectMap>& cfp = map_.cellsFromPointsMap();

        forAll(cfp, cfpI)
        {
            cm[cfp[cfpI].index()] = 0;
        }

        const List<objectMap>& cfe = map_.cellsFromEdgesMap();

        forAll(cfe, cfeI)
        {
            cm[cfe[cfeI].index()] = 0;
        }

        const List<objectMap>& cff = map_.cellsFromFacesMap();

        forAll(cff, cffI)
        {
            cm[cff[cffI].index()] = 0;
        }

        const List<objectMap>& cfc = map_.cellsFromCellsMap();

        forAll(cfc, cfcI)
        {
            cm[cfc[cfcI].index()] = 0;
        }

        if (min(cm) < 0)
        {
            insertedCells_ = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellMapper::~cellMapper()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellMapper::sizeBeforeMapping() const
{
    return map_.nOldCells();
}


const Foam::labelUList& Foam::cellMapper::directAddressing() const
{
    if (!direct())
    {
        FatalErrorInFunction
            << "Requested direct addressing for an interpolative mapper."
            << abort(FatalError);
    }

    if (!insertedObjects())
    {
        // No inserted cells.  Re-use cellMap
        return map_.cellMap();
    }
    else
    {
        if (!directAddrPtr_)
        {
            calcAddressing();
        }

        return *directAddrPtr_;
    }
}


const Foam::labelListList& Foam::cellMapper::addressing() const
{
    if (direct())
    {
        FatalErrorInFunction
            << "Requested interpolative addressing for a direct mapper."
            << abort(FatalError);
    }

    if (!interpolationAddrPtr_)
    {
        calcAddressing();
    }

    return *interpolationAddrPtr_;
}


const Foam::scalarListList& Foam::cellMapper::weights() const
{
    if (direct())
    {
        FatalErrorInFunction
            << "Requested interpolative weights for a direct mapper."
            << abort(FatalError);
    }

    if (!weightsPtr_)
    {
        calcAddressing();
    }

    return *weightsPtr_;
}


const Foam::labelList& Foam::cellMapper::insertedObjectLabels() const
{
    if (!insertedCellLabelsPtr_)
    {
        if (!insertedObjects())
        {
            // There are no inserted cells
            insertedCellLabelsPtr_ = new labelList(0);
        }
        else
        {
            calcAddressing();
        }
    }

    return *insertedCellLabelsPtr_;
}


// ************************************************************************* //
