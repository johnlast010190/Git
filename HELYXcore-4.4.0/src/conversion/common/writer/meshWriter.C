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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "common/writer/meshWriter.H"
#include "meshes/meshShapes/cellModeller/cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::meshWriter::defaultMeshName = "meshExport";


const Foam::cellModel* Foam::meshWriter::unknownModel = Foam::cellModeller::
lookup
(
    "unknown"
);


const Foam::cellModel* Foam::meshWriter::tetModel = Foam::cellModeller::
lookup
(
    "tet"
);


const Foam::cellModel* Foam::meshWriter::pyrModel = Foam::cellModeller::
lookup
(
    "pyr"
);


const Foam::cellModel* Foam::meshWriter::prismModel = Foam::cellModeller::
lookup
(
    "prism"
);


const Foam::cellModel* Foam::meshWriter::hexModel = Foam::cellModeller::
lookup
(
    "hex"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriter::meshWriter
(
    const polyMesh& mesh,
    const scalar scaling
)
:
    mesh_(mesh),
    scaleFactor_(scaling),
    boundaryRegion_(),
    cellTable_(),
    cellTableId_()
{
    // Sanity
    if (scaleFactor_ <= VSMALL)
    {
        scaleFactor_ = 1;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriter::~meshWriter()
{}


// ************************************************************************* //
