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
    Foam::functionObjects::runTimeVis::PvdReaderDataSetProvider

Description
    Provider for the objects that are saved as pvd files

SourceFiles
    PvdReaderDataSetProvider.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <utility>

#include "itemDataSetProviders/itemDataSetProvider.H"
#include "dataStructs/fileSource/pvdFileData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class engysPVDReader;

namespace Foam::functionObjects::runTimeVis
{


/*---------------------------------------------------------------------------*\
              Class PvdReaderDataSetProvider Declaration
\*---------------------------------------------------------------------------*/

class PvdReaderDataSetProvider : public ItemDataSetProvider
{
private:
    PvdFileData data_;
    vtkSmartPointer<engysPVDReader> reader_;

public:
    PvdReaderDataSetProvider(const std::string &name, PvdFileData data);
    ~PvdReaderDataSetProvider() override = default;

protected:
    void update(scalar currentTime) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
