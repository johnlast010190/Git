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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    faPatchFieldMapperPatchRef

Description

\*---------------------------------------------------------------------------*/

#ifndef faPatchFieldMapperPatchRef_H
#define faPatchFieldMapperPatchRef_H

#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class faPatchFieldMapperPatchRef Declaration
\*---------------------------------------------------------------------------*/

class faPatchFieldMapperPatchRef
:
    public faPatchFieldMapper
{
    // Private data

        const faPatch& sourcePatch_;
        const faPatch& targetPatch_;


public:

    // Constructors

        //- Construct from components
        faPatchFieldMapperPatchRef(const faPatch& source, const faPatch& target)
        :
            sourcePatch_(source),
            targetPatch_(target)
        {}

        //- Disallow default bitwise copy construction
        faPatchFieldMapperPatchRef(const faPatchFieldMapperPatchRef&) = delete;


    // Member functions

        const faPatch& sourcePatch() const
        {
            return sourcePatch_;
        }

        const faPatch& targetPatch() const
        {
            return targetPatch_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faPatchFieldMapperPatchRef&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
