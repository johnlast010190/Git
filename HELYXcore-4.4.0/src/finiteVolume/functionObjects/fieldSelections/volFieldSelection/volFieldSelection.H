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
    (c) 2017-19 OpenCFD Ltd.

Class
    Foam::functionObjects::volFieldSelection

Description
    Helper class to manage volume field selections

SourceFiles
    volFieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volFieldSelection_H
#define functionObjects_volFieldSelection_H

#include "functionObjects/fieldSelections/fieldSelection/fieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class volFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class volFieldSelection
:
    public fieldSelection
{
private:

    // Private Member Functions

        //- No copy construct
        volFieldSelection(const volFieldSelection&) = delete;


protected:

    // Protected Member Functions

        //- Add registered GeometricField types to selection
        template<template<class> class PatchType, class MeshType>
        void addRegisteredGeoFields(DynamicList<fieldInfo>& set) const;


public:

    //- Construct from object registry
    volFieldSelection
    (
        const objectRegistry& obr,
        const bool includeComponents = false
    );


    //- Destructor
    virtual ~volFieldSelection() = default;


    // Member Functions

        //- Update the selection using current contents of obr_
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionObjects/fieldSelections/volFieldSelection/volFieldSelectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
