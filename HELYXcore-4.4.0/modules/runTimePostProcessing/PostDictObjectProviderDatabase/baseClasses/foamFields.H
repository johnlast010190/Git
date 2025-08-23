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
    (c) 2019 OpenCFD Ltd.
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::FoamFields

Description
    Class that handles sets of foamFields

SourceFiles
    <FoamField.C>

\*---------------------------------------------------------------------------*/
#pragma once

#include "primitives/strings/stringOps/stringOps.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "foamField.H"
#include <unordered_set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class FoamFields
{
public:

    // Constructors
    FoamFields() = default;
    ~FoamFields() = default;

    void addField(const foamField& field);

    void addFoamField(const foamField& field);

    void addAllFieldsMarker();

    /**
     * @return a vector with the Foam fields, which are just the field names (eg. U, p)
     */
    [[nodiscard]] std::vector<foamField> getFoamFields() const;

    /**
     * @return returns a list of all foam fields that require cell data but not point data
     */
    [[nodiscard]] std::vector<foamField> getCellOnlyFoamFields() const;

    /**
     * @return a vector with the color fields, which are the fields and the component markers (eg. U-X, p)
     */
    [[nodiscard]] std::vector<foamField> getColorFields() const;

    /**
     * @return returns true if at least one field as a point association
     */
    [[nodiscard]] bool hasPointField() const;

    /**
     * @return a set with the foam fields that have a point/cell association
     */
    [[nodiscard]] std::set<std::string> getPointFoamFieldsSet() const;
    [[nodiscard]] std::set<std::string> getCellFoamFieldsSet() const;

    void merge(const FoamFields& foamFields);
    void substituteAllFieldsMarker(const FoamFields& allFields);
    void clear();

    [[nodiscard]] bool hasAllFieldsMarker() const;
    [[nodiscard]] bool containsFoamField(const foamField& foam) const;
    [[nodiscard]] bool empty() const;

private:

    struct FoamFieldHash
    {
        std::size_t operator () (const foamField& field) const
        {
            return std::hash<std::string>()(field);
        }
    };

    typedef std::unordered_set<foamField, FoamFieldHash> FieldsSet;

    FieldsSet fieldsSet;
};

} // End namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
