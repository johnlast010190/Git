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
    (c) 2014 Engys Ltd.

Class
    Foam::modifiedNewtonsMethod

Description
    Class for calculating the optimization step given the Hessian tensor and the
    derivative. The method performs a modified Cholesky decomposition and a tensor
    inversion

SourceFiles
      modifiedNewtonsMethod.C
\*---------------------------------------------------------------------------*/

#ifndef modifiedNewtonsMethod_H
#define modifiedNewtonsMethod_H

#include "primitives/Tensor/tensor/tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class modifiedNewtonsMethod Declaration
\*---------------------------------------------------------------------------*/

class modifiedNewtonsMethod
{

public:

    // Constructors
        modifiedNewtonsMethod();

        //- Destructor
        ~modifiedNewtonsMethod()
        {}
    // Member Functions

        static void modifiedCholesky(tensor& T);

        static vector newtonsStep(const vector& G, const tensor& T);
};

} //End namespace Foam

#endif
