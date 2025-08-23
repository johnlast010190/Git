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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fv::cellPartialLimitedGrad

Group
    grpFvGradSchemes

Description

SourceFiles
    cellPartialLimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef cellPartialLimitedGrad_H
#define cellPartialLimitedGrad_H

#include "finiteVolume/gradSchemes/gradScheme/gradScheme.H"
#include "finiteVolume/gradSchemes/limitedGradSchemes/cellLimitedGrad/cellLimitedGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class cellPartialLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellPartialLimitedGrad
:
    public fv::cellLimitedGrad<Type>
{
    // Private Data

        //- Name of meshObject
        const word meshObjectName_;

    // Private Functions
        virtual void blendLimiter(Field<Type>&) const;

        //- write the cells that will be limited - only for debugging
        void writeField(const boolList&) const;

public:

    //- RunTime type information
    TypeName("cellPartialLimited");


    // Constructors

        //- Construct from mesh and schemeData
        cellPartialLimitedGrad
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            cellLimitedGrad<Type>(mesh, db, schemeData),
            meshObjectName_(schemeData)
        {
        }

        //- Disallow default bitwise copy construction
        cellPartialLimitedGrad(const cellPartialLimitedGrad&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellPartialLimitedGrad&) = delete;
};


// * * * * * * * * * * * * Inline Member Function  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteVolume/gradSchemes/limitedGradSchemes/cellPartialLimitedGrad/cellPartialLimitedGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
