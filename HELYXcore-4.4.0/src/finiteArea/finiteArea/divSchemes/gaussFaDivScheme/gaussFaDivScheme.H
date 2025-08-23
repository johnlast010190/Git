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
    (c) 2016-2024 Engys Ltd.

Class
    gaussDivScheme

Description
    Basic second-order div using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaDivScheme.C
    gaussFaDivSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gaussFaDivScheme_H
#define gaussFaDivScheme_H

#include "finiteArea/divSchemes/faDivScheme/faDivScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class gaussDivScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussDivScheme
:
    public fa::divScheme<Type>
{
public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussDivScheme(const faMesh& mesh)
        :
            divScheme<Type>(mesh)
        {}

        //- Construct from Istream
        gaussDivScheme(const faMesh& mesh, Istream& is)
        :
            divScheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        gaussDivScheme(const gaussDivScheme&) = delete;


    //- Destructor
    virtual ~gaussDivScheme()
    {}


    // Member Functions

        tmp
        <
            GeometricField
            <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
        > facDiv
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussDivScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/divSchemes/gaussFaDivScheme/gaussFaDivScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
