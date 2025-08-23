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
    Foam::functionObjects::runTimeVis::ExternalFields

Description
    contains references to the foam meshes, and functions to access their properties

SourceFiles
    externalFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimeVis_externalFields
#define functionObjects_runTimeVis_externalFields

#include "itemDataSetProviders/externalItemDataSetProvider.H"

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "primitives/ranges/MinMax/MinMax.H"

#include "baseClasses/foamFields.H"
#include "baseClasses/foamField.H"
#include "dictionaries/dictionaries.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

/*---------------------------------------------------------------------------*\
                        Class ExternalFields Declaration
\*---------------------------------------------------------------------------*/
class ExternalFields
{
private:
    std::vector<ExternalItemDataSetProvider*> providers_;

    HashTable<scalarMinMax, word> domainRanges_;
    FoamFields colourFields_;

public:
    ExternalFields();
    ~ExternalFields() = default;

    void addColourField(const foamField& colourField);

    void updateForTimeStep(label timeIndex);

    [[nodiscard]] scalarMinMax getDomainRangeForField(const foamField& field) const;

    void insertProvider(ExternalItemDataSetProvider *p);
    void clearProviders();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
