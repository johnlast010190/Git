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
    (c) 2020-2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::foamField

Description
    Class that deals intelligently with Helyx GUI field naming conventions.

SourceFiles
    <foamField.C>

\*---------------------------------------------------------------------------*/
#pragma once

#include "primitives/strings/stringOps/stringOps.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include <unordered_map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects::runTimeVis
{

class foamField : public string
{
protected:
    static const char *INDEX_COLOR;
    static const char *SOLID_COLOR;
    static const char *SOLID_ASSOCIATION;
    static const char *INDEXED_ASSOCIATION;
    static const char *POINT_ASSOCIATION;
    static const char *CELL_ASSOCIATION;
    static const char *MAG_COMPONENT;
    static const char *ALL_FIELDS_KEY;

    static const std::set<std::string>& getAllowedAssociations();
    static const std::unordered_map<std::string, int>& getAllowedComponents();

    // The name as it appears in OpenFOAM
    string foamName_ = SOLID_COLOR;

    // The name without field association, but with component specifier (equal
    // to foamName_ for scalars)
    string nameLessAssociation_ = SOLID_COLOR;

    word component_ = "";
    word association_ = SOLID_ASSOCIATION;
    int componentIndex_ = 0;

    bool isComponent_ = false;
    bool isMagnitude_ = false;
    bool isSolidColour_ = false;

    [[nodiscard]] string lessExt() const;

    void initialise();


public:

    // Constructors
    foamField() = default;
    ~foamField() = default;

    // Construct from string
    explicit foamField(const string &fullName) : string(fullName)
    {
        initialise();
    };

    explicit foamField(const std::string &fullName) : string(fullName)
    {
        initialise();
    };

    static foamField withCellAssociation(const std::string &foamName);
    static foamField withPointAssociation(const std::string &foamName);

    [[nodiscard]] static foamField AllFieldsMarker()
    {
        return foamField(std::string(ALL_FIELDS_KEY));
    }

    // Accessors
    [[nodiscard]] inline string getFoamName() const { return foamName_; }
    [[nodiscard]] inline word getComponent() const { return component_; }
    [[nodiscard]] inline label getComponentIndex() const { return componentIndex_; }
    [[nodiscard]] inline word getAssociation() const { return association_; }
    [[nodiscard]] inline word lessAssociation() const { return nameLessAssociation_; }
    [[nodiscard]] inline bool isComponent() const { return isComponent_; }
    [[nodiscard]] inline bool isSolidColour() const { return isSolidColour_; }
    [[nodiscard]] inline bool isIndexedColour() const { return isSolidColour_ && INDEXED_ASSOCIATION == association_; }
    [[nodiscard]] inline bool isPointAssociation() const { return POINT_ASSOCIATION == association_; }
    [[nodiscard]] inline bool isCellAssociation() const { return CELL_ASSOCIATION == association_; }
    [[nodiscard]] inline bool hasFieldAssociation() const { return isPointAssociation() || isCellAssociation(); }
    [[nodiscard]] inline bool isMagnitude() const { return isMagnitude_; }
    [[nodiscard]] inline bool isAllFieldsMarker() const { return *this == ALL_FIELDS_KEY; }


    friend Ostream &operator<<(Ostream &os, const foamField &fn)
    {
        os << string(fn.data())
           << " (Component: " << fn.component_
           << ", foam field name: " << fn.foamName_
           << ", association: \"" << fn.association_ << "\""
           << ")";
        return os;
    };

    friend Istream &operator>>(Istream &is, foamField &m)
    {
        string w;
        is >> w;
        m = foamField(w);
        return is;
    }
};

} // End namespace Foam::functionObjects::runTimeVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
