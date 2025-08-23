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
    (c) 2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "extractFromFile.H"
#include "edgeMesh/edgeMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{
    addNamedToRunTimeSelectionTable
    (
        method,
        extractFromFile,
        dictionary,
        extractFromFile
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::extractFromFile::extractFromFile
(
    const dictionary& dict
)
:
    method()
{
    const dictionary& coeffDict =
        dict.optionalSubDict("extractFromFileCoeffs");

    featureEdgeFile_ = coeffDict.lookup<fileName>("featureEdgeFile");
    coeffDict.readIfPresent("geometricTestOnly", geometricTestOnly_);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::extractFromFile::~extractFromFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::autoPtr<Foam::surfaceFeatures>
Foam::surfaceFeaturesExtraction::extractFromFile::features
(
    const triSurface& surf
) const
{
    edgeMesh eMesh(featureEdgeFile_);

    // Sometimes duplicate edges are present. Remove them.
    eMesh.mergeEdges();

    Info<< nl << "Reading existing feature edges from file "
        << featureEdgeFile_ << nl
        << "Selecting edges based purely on geometric tests: "
        << geometricTestOnly().asText() << endl;

    return autoPtr<surfaceFeatures>
    (
        new surfaceFeatures
        (
            surf,
            eMesh.points(),
            eMesh.edges(),
            1e-6,  // mergeTol
            geometricTestOnly()
        )
    );
}


// ************************************************************************* //
