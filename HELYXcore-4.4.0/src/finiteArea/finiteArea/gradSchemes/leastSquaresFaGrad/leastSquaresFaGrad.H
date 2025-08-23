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
    (c) 2017 Wikki Ltd.
    (c) 2016-2024 Engys Ltd.

Class
    Foam::fa::leastSquaresFaGrad

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    leastSquaresFaGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaGrad_H
#define leastSquaresFaGrad_H

#include "finiteArea/gradSchemes/faGradScheme/faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresFaGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresFaGrad
:
    public fa::gradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("leastSquares");


    // Constructors

        //- Construct from mesh
        leastSquaresFaGrad(const faMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresFaGrad(const faMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}

        //- Disallow default bitwise copy construction
        leastSquaresFaGrad(const leastSquaresFaGrad&) = delete;


    //- Destructor
    virtual ~leastSquaresFaGrad()
    {}


    // Member Functions

        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        >grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresFaGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/gradSchemes/leastSquaresFaGrad/leastSquaresFaGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
