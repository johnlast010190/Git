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
    (c) 2019 OpenCFD Ltd.
    (c) 2022 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::Font

Description
    Enum with the possible fonts for the color legend

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "containers/HashTables/HashTable/HashTable.H"

namespace Foam::functionObjects::runTimeVis
{
    class ItemType
    {
    public:
        enum Value
        {
            PATCH,
            BOUNDARY_GROUP = PATCH, // externalBoundaries is considered a patch
            SURFACE,
            GROUP,
            SLICE,
            CLIP,
            TRANSFORM,
            ISO_SURFACE,
            STREAMLINES,
            GLYPH,
            THRESHOLD,
            FIELD_SAMPLING,
            FIELD_SAMPLING_FILE_SOURCE,
            CELLZONE,
            FACEZONE,
            PROCESSOR_BOUNDARY,
            FEATURE_LINE,
            LINE_3D,
            MERIDIONAL_GRID,
            TURBO_SLICE_STREAMWISE,
            TURBO_SLICE_SPANWISE,
            INLET_TO_OUTLET_LINE,
            HUB_TO_SHROUD_LINE,
            BLADE_LOADING_LINE,
            PVD_FILE_SOURCE,
            IMPORT_DATA,
            TOOL_3D_TO_2D,
            SURFACE_SPLITTER_TEST
        };

        ItemType(){ value = PATCH; }
        explicit ItemType(Value aValue) {value = aValue;}
        explicit ItemType(const word& aType) { value = dictItemTypes_[aType]; }
        ItemType(const ItemType& type) {value = type.value;}

        [[nodiscard]] int getValue() const {return value;}

        [[nodiscard]] word getName() const
        {
            if (value == PROCESSOR_BOUNDARY)
            {
                return "processorBoundary";
            }
            for (auto itemType = dictItemTypes_.begin(); itemType != dictItemTypes_.end(); itemType++)
            {
                if (itemType.object() == value)
                {
                    return itemType.key();
                }
            }
            return "unknown";
        }

        // Due to a bug in GCC 5.1, these cannot be constexpr
        bool operator==(const ItemType& i) const { return value == i.value; }
        bool operator!=(const ItemType& i) const { return value != i.value; }

        ItemType& operator=(const ItemType& i) {value = i.value; return *this;}
        ItemType& operator=(const Value& i) {value = i; return *this;}

        friend Ostream& operator<<(Ostream& os, const ItemType& itemType)
        {
            os << itemType.getName();
            return os;
        };

    private:
        HashTable<Value, word> dictItemTypes_
        {
            {"patch"            , PATCH             },
            {"surface"          , SURFACE           },
            {"group"            , GROUP             },
            {"slice"            , SLICE             },
            {"clip"             , CLIP              },
            {"transform"        , TRANSFORM         },
            {"isoSurface"       , ISO_SURFACE       },
            {"streamlines"      , STREAMLINES       },
            {"glyph"            ,       GLYPH             },
            {"threshold"        ,       THRESHOLD         },
            {"fieldSampling"    ,       FIELD_SAMPLING    },
            {"fieldSamplingFileSource", FIELD_SAMPLING_FILE_SOURCE},
            {"cellZone"         ,       CELLZONE          },
            {"faceZone"         ,       FACEZONE          },
            {"featureline"      ,       FEATURE_LINE      },
            {"line3d"           ,       LINE_3D           },
            {"turboMeridionalGrid",     MERIDIONAL_GRID   },
            {"turboSliceStreamwise",    TURBO_SLICE_STREAMWISE},
            {"turboSliceSpanwise",      TURBO_SLICE_SPANWISE},
            {"bladeLoadingLine" ,       BLADE_LOADING_LINE},
            {"inletToOutletLine",       INLET_TO_OUTLET_LINE},
            {"hubToShroudLine"  ,       HUB_TO_SHROUD_LINE},
            {"pvdFileSource"    ,       PVD_FILE_SOURCE},
            {"importData"       ,       IMPORT_DATA},
            {"tool3Dto2D"       ,       TOOL_3D_TO_2D},
            {"surfaceSplitterTest",     SURFACE_SPLITTER_TEST}
        };
        Value value;
    } ;
} // End namespace Foam
