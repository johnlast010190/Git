/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.0.1
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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::IsoSurfaceDataSetProvider

Description
    Provider for the iso surface 3d objects, which computes their vtk dataset

SourceFiles
    isosurfaceDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "itemDataSetProviders/itemDataSetProvider.H"

#include "hash/hasher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkContourFilter;

namespace Foam::functionObjects::runTimeVis
{

struct IsoSurfaceObjectData
{
    List<scalar> values;
    foamField scalarField;
    bool enableSmoothing;
    label iterations;
    scalar convergence;

    bool operator==(const IsoSurfaceObjectData& other) const
    {
        if (values.size() != other.values.size())
        {
            return false;
        }
        bool res = scalarField == other.scalarField &&
               enableSmoothing == other.enableSmoothing &&
                    iterations == other.iterations &&
                   convergence == other.convergence;
        if (!res)
        {
            return false;
        }
        return std::all_of(values.begin(), values.end(),
                           [&other](scalar v) {
                                return other.values.found(v);
                            }
                        );
    }

    [[nodiscard]] size_t computeAndReturnHash() const
    {
        size_t hash = 0;
        hasher::hash_combine(hash, scalarField);
        hasher::hash_combine(hash, enableSmoothing);
        hasher::hash_combine(hash, iterations);
        hasher::hash_combine(hash, convergence);

        for (scalar v : values) {
            hasher::hash_add_unordered(hash, v);
        }
        return hash;
    }
};

/*---------------------------------------------------------------------------*\
                 Class TransformDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class IsoSurfaceDataSetProvider : public ItemDataSetProvider
{
private:

protected:

        vtkSmartPointer<vtkContourFilter> contourFilter_;

        void update(scalar currentTime) override;

public:
        IsoSurfaceDataSetProvider(const std::string& name, const IsoSurfaceObjectData& dictData);
        ~IsoSurfaceDataSetProvider() override = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
