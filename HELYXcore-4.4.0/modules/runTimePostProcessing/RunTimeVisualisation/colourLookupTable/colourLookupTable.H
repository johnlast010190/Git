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

Description
    Class that generates VTK colour lookup tables for runTimeVisualisation
    scenes.

Class
    Foam::functionObjects::ColourLookupTable

\*---------------------------------------------------------------------------*/

#pragma once

#include "colourMaps/colourMaps.H"
#include "dataStructs/scene/colourLookupTableData.H"
#include "baseClasses/foamField.H"

#include "vtkLookupTable.h"
#include "vtkDiscretizableColorTransferFunction.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                   Class ColourLookupTable Declaration
\*---------------------------------------------------------------------------*/

// Takes input from volMeshFilter and generates a vtkLookupTable of colours
class ColourLookupTable
{
protected:

    // Protected data

        ColourLookupTableData data_;
        scalarMinMax previousDomainMinMax_;

        vtkSmartPointer<vtkLookupTable> actorLookupTable_;
        vtkSmartPointer<vtkDiscretizableColorTransferFunction> legendLookupTable_;

        vtkSmartPointer<vtkDiscretizableColorTransferFunction> colorTransferFunction_;

    // Protected member functions

        void initializeColorTransferFunction
        (
                const ColourMaps& colourMaps
        );

        void initializeActorLookupTable();
        void initializeLegendLookupTable();

        void setSharedLookupTableProperties
        (
                vtkScalarsToColors* lookupTable
        ) const;

        static scalarMinMax adjustRange(scalarMinMax range);

        void buildActorLookupTable(scalarMinMax range);
        void buildLegendLookupTable(scalarMinMax range);

        bool mustUpdateRange(scalarMinMax domainMinMax);


public:

    // Constructors

        ColourLookupTable
        (
            const ColourLookupTableData& data,
            const ColourMaps& colourMaps
        );


    // Member Functions

        void updateDomainRange(scalarMinMax domainMinMax);

        // Access

        [[nodiscard]] vtkLookupTable* getActorLookupTable() const
        {
            return actorLookupTable_;
        }

        [[nodiscard]] vtkDiscretizableColorTransferFunction* getLegendLookupTable() const
        {
            return legendLookupTable_;
        }

    void setNanTableProperties();
    void setColorTableProperties(scalarMinMax &adjustedRange);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace


// ************************************************************************* //
