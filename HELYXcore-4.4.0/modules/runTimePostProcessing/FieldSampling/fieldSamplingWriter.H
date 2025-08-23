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
    Foam::functionObjects::fieldSampling::FieldSamplingWriter

Description
    Resamples the volume mesh field values to a uniform grid

SourceFiles
    fieldSamplingWriter.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/fileName/fileName.H"
#include "vtkSmartPointer.h"
#include "db/Time/Time.H"
#include "db/functionObjects/writeFile/writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysResampleToImage;
class engysPVDWriter;
class vtkDataSetAttributes;

namespace Foam::functionObjects
{

namespace runTimeVis
{
class PostDictObjectProviderDatabase;
}

namespace fieldSample
{

struct FieldSamplingData;

/*---------------------------------------------------------------------------*\
                  Class FieldSamplingWriter Declaration
\*---------------------------------------------------------------------------*/

class FieldSamplingWriter
{
private:

    vtkSmartPointer<engysResampleToImage> fieldSamplingFilter_;
    vtkSmartPointer<engysPVDWriter> pvdWriter_;

    const FieldSamplingData &data_;

    const Switch& log;

    void initialiseSampler();

public:
    explicit FieldSamplingWriter(const FieldSamplingData &data, const Switch& log);

    void write(
        const Time &currentTime,
        const runTimeVis::PostDictObjectProviderDatabase &database
    );

    void clearSampler();

    ~FieldSamplingWriter() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldSample
} // End namespace Foam::functionObjects

// ************************************************************************* //
