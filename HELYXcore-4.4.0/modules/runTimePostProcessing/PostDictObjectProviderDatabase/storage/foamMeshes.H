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
    (c) 2022-2024 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FoamMeshes

Description
    contains references to the foam meshes, and functions to access their properties

SourceFiles
    foamMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_foamMeshes_H
#define functionObjects_runTimeVis_foamMeshes_H

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "primitives/ranges/MinMax/MinMax.H"

#include "baseClasses/foamFields.H"
#include "baseClasses/foamField.H"
#include "baseClasses/meshAndFields.H"
#include "dictionaries/dictionaries.H"

#include <unordered_set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                        Class ItemStorage Declaration
\*---------------------------------------------------------------------------*/
class FoamMeshes
{
private:
    HashTable<MeshAndFields, word> meshes_;

    HashTable<scalarMinMax, word> domainRanges_;
    boundBox domainBounds_;

    FoamFields colourFields_;

public:
    //- Gets the foam meshes from the core
    FoamMeshes(const Time& time, const Dictionaries &dictionaries);

    void addColourField(const foamField& colourField);

    void updateDomainRangeForColourFields();
    void updateDomainBounds();

    [[nodiscard]] scalarMinMax getDomainRangeForField(const foamField& field) const;

    [[nodiscard]] boundBox getDomainBounds() const {return domainBounds_;};

    [[nodiscard]] MeshAndFields getMesh(const word& region) const {return meshes_[region];};
    [[nodiscard]] MeshAndFields getDefaultMesh() const {return *meshes_.begin();};
    [[nodiscard]] const HashTable<MeshAndFields, word>& getMeshes() const {return meshes_;};

private:
    void addFoamMesh(
        const Time &time,
        const word &meshRegionName,
        const word& instanceName
    );

    static const objectRegistry &getFieldObjectRegistry(
        const Time &time,
        const word &meshRegionName,
        const word &instanceName
    );
    void updateDomainRangeForList(const FoamFields &fields);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
