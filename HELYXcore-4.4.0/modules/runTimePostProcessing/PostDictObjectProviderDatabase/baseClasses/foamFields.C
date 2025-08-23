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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "baseClasses/foamFields.H"

namespace Foam::functionObjects::runTimeVis
{

// ************************************************************************* //
void FoamFields::addField(const foamField &field)
{
    fieldsSet.insert(field);
}

void FoamFields::addFoamField(const foamField &field)
{
    if (field.hasFieldAssociation() || field.isAllFieldsMarker())
    {
        fieldsSet.insert(field);
    }
    else
    {
        fieldsSet.insert(foamField::withPointAssociation(field));
        fieldsSet.insert(foamField::withCellAssociation(field));
    }
}

//----------------------------------------------------------------------
void FoamFields::addAllFieldsMarker()
{
    addField(foamField::AllFieldsMarker());
}

bool compFoamFields(const foamField& v1, const foamField& v2)
{
    size_t i = 0;
    const size_t size1 = v1.size();
    const size_t size2 = v2.size();
    while (i < size1 && i < size2)
    {
        unsigned char c1 = std::tolower(v1[i]);
        unsigned char c2 = std::tolower(v2[i]);
        if (c1 != c2) return c1 < c2;
        i++;
    }
    return size1 < size2;
}

//----------------------------------------------------------------------
std::vector<foamField> FoamFields::getFoamFields() const
{
    FieldsSet foamFieldsSet;
    for (const foamField& field : fieldsSet)
    {
        if (field.hasFieldAssociation())
        {
            foamFieldsSet.insert(foamField(field.getFoamName()));
        }
    }
    std::vector<foamField> sortedFields(foamFieldsSet.begin(), foamFieldsSet.end());
    std::sort(sortedFields.begin(), sortedFields.end(), compFoamFields);
    return sortedFields;
}

//----------------------------------------------------------------------
std::vector<foamField> FoamFields::getCellOnlyFoamFields() const
{
    FieldsSet foamFieldsSet;
    for (const foamField& field : fieldsSet)
    {
        if (field.isCellAssociation())
        {
            foamFieldsSet.insert(foamField(field.getFoamName()));
        }
    }
    for (const foamField& field : fieldsSet)
    {
        if (field.isPointAssociation())
        {
            foamFieldsSet.erase(foamField(field.getFoamName()));
        }
    }
    std::vector<foamField> cellFields(foamFieldsSet.begin(), foamFieldsSet.end());
    return cellFields;
}

//----------------------------------------------------------------------
std::vector<foamField> FoamFields::getColorFields() const
{
    FieldsSet foamFieldsSet;
    for (const foamField& field : fieldsSet)
    {
        foamFieldsSet.insert(foamField(field.lessAssociation()));
    }
    return {foamFieldsSet.begin(), foamFieldsSet.end()};
}

//----------------------------------------------------------------------
bool FoamFields::hasPointField() const
{
    return std::any_of
    (
        fieldsSet.begin(),
        fieldsSet.end(),
        [](const foamField &field) { return field.isPointAssociation(); }
    );
}

//----------------------------------------------------------------------
bool FoamFields::hasAllFieldsMarker() const
{
    return fieldsSet.find(foamField::AllFieldsMarker()) != fieldsSet.end();
}

//----------------------------------------------------------------------
void FoamFields::merge(const FoamFields& foamFields)
{
    fieldsSet.insert(foamFields.fieldsSet.begin(), foamFields.fieldsSet.end());
}

//----------------------------------------------------------------------
void FoamFields::substituteAllFieldsMarker(const FoamFields& allFields)
{
    merge(allFields);
    fieldsSet.erase(foamField::AllFieldsMarker());
}

//----------------------------------------------------------------------
void FoamFields::clear()
{
    fieldsSet.clear();
}

//----------------------------------------------------------------------
bool FoamFields::containsFoamField(const foamField &foam) const
{
    return std::any_of(fieldsSet.begin(), fieldsSet.end(),
                       [foam](const foamField& field) -> bool {return field.getFoamName() == foam.getFoamName();}
                       );
}

//----------------------------------------------------------------------
bool FoamFields::empty() const
{
    return fieldsSet.empty();
}

//----------------------------------------------------------------------
std::set<std::string> FoamFields::getPointFoamFieldsSet() const
{
    std::set<std::string> foamFieldsSet;
    for (const foamField& field : fieldsSet)
    {
        if (field.isPointAssociation()) foamFieldsSet.insert(field.getFoamName());
    }
    return foamFieldsSet;
}

//----------------------------------------------------------------------
std::set<std::string> FoamFields::getCellFoamFieldsSet() const
{
    std::set<std::string> foamFieldsSet;
    for (const foamField& field : fieldsSet)
    {
        if (field.isCellAssociation()) foamFieldsSet.insert(field.getFoamName());
    }
    return foamFieldsSet;
}

} // namespace
