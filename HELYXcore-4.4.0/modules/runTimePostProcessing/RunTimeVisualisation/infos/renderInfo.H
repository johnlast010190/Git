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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::RenderInfo

Description
    A struct with basic info about the rendering

SourceFiles
    <none>

\*---------------------------------------------------------------------------*/

#pragma once

#include "baseClasses/foamField.H"
#include "types/pvdFieldsType.H"

#include "primitives/ints/label/label.H"
#include "containers/Lists/List/List.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

enum EDFCompressionLevel
{
    FASTEST,
    FAST,
    BALANCED,
    SMALL,
    SMALLEST
};

struct RenderInfo {
    bool exportPng = false;
    bool exportX3d = false;
    bool exportEdf = false;
    label width = 0;
    label height = 0;
    bool transparentBackground = false;
    bool cropToContents = false;
    bool centerToCroppedContents = false;
    List<word> pvdFieldsList;
    PvdFieldsType pvdFieldsType;
    List<foamField> edfFields;
    EDFCompressionLevel edfCompressionLevel;

    explicit RenderInfo(const dictionary& postDict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
