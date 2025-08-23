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
    (c) 2015-2019 OpenCFD Ltd.
    (c) 2020-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "baseClasses/foamField.H"

namespace Foam::functionObjects::runTimeVis
{
const char *foamField::INDEX_COLOR = "index";
const char *foamField::SOLID_COLOR = "solidColor";
const char *foamField::SOLID_ASSOCIATION = "SOLID";
const char *foamField::INDEXED_ASSOCIATION = "INDEXED";
const char *foamField::POINT_ASSOCIATION = "POINT";
const char *foamField::CELL_ASSOCIATION = "CELL";
const char *foamField::MAG_COMPONENT = "Mag";
const char *foamField::ALL_FIELDS_KEY = "ALL_FIELDS_MARKER";

const std::set<std::string>& Foam::functionObjects::runTimeVis::foamField::getAllowedAssociations()
{
    const static std::set<std::string> ALLOWED_ASSOCIATIONS = {POINT_ASSOCIATION,
                                                               CELL_ASSOCIATION,
                                                               SOLID_ASSOCIATION,
                                                               INDEXED_ASSOCIATION};
    return ALLOWED_ASSOCIATIONS;
}

const std::unordered_map<std::string, int>& Foam::functionObjects::runTimeVis::foamField::getAllowedComponents()
{
    const static std::unordered_map<std::string, int> ALLOWED_COMPONENTS = {
        {MAG_COMPONENT, -1},
        {"X", 0},
        {"Y", 1},
        {"Z", 2},
        {"XX", 0},
        {"YY", 1},
        {"ZZ", 2},
        {"XY", 3},
        {"XZ", 4},
        {"YZ", 5}
    };
    return ALLOWED_COMPONENTS;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

string foamField::lessExt() const
{
    const size_type i = find_ext();
    return substr(0, i);
}

void foamField::initialise()
{
    association_ = SOLID_ASSOCIATION;
    nameLessAssociation_ = this->c_str();
    if (hasExt())
    {
        word extension = ext();

        const std::set<std::string>& allowedAssociations = getAllowedAssociations();
        if (allowedAssociations.find(extension) != allowedAssociations.end())
        {
            association_ = extension;
            nameLessAssociation_ = this->lessExt();
        }
    }

    SubStrings<string> splitFieldName = stringOps::split(nameLessAssociation_, '-');
    if (splitFieldName.size() > 1)
    {
        word potentialComponent(splitFieldName.last());
        const auto& allowedComponents = getAllowedComponents();
        const auto& itr = allowedComponents.find(potentialComponent);
        if (itr != allowedComponents.end())
        {
            // Handle vector component specifier
            component_ = potentialComponent;
            foamName_ = substr(0, find_last_of("-"));
            componentIndex_ = itr->second;
            if (componentIndex_ < 0)
            {
                isMagnitude_ = true;
            }
            else
            {
                isComponent_ = true;
            }
        }
        else
        {
            // it's possible to fall here if the field has a name like alpha.air-ISO
            foamName_ = nameLessAssociation_;
        }
    }
    else  // A scalar field
    {
        foamName_ = nameLessAssociation_;
    }

    if (SOLID_COLOR == foamName_ || INDEX_COLOR == foamName_)
    {
        isSolidColour_ = true;
        const word expectedAssociation = SOLID_COLOR == foamName_ ? SOLID_ASSOCIATION : INDEXED_ASSOCIATION;
        if (expectedAssociation != association_)
        {
            WarningInFunction
                << "Field with foam name \""
                << foamName_ << "\" had unexpected association \""
                << association_
                << "\" (was expecting " << expectedAssociation << "). "
                << "Defaulting to " << expectedAssociation << "."
                << endl;
        }
        association_ = expectedAssociation;
    }
}

foamField foamField::withPointAssociation(const std::string& foamName)
{
    foamField field(foamName);
    if (field.isPointAssociation()) return field;
    return foamField{field.lessAssociation() + "." + POINT_ASSOCIATION};
}

foamField foamField::withCellAssociation(const std::string& foamName)
{
    foamField field(foamName);
    if (field.isCellAssociation()) return field;
    return foamField{field.lessAssociation() + "." + CELL_ASSOCIATION};
}

// ************************************************************************* //
} // namespace