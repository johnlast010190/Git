/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Class
    Foam::forceEquation

Description
    Force a fvMatrix to fixed values in specific places

SourceFiles
    forceEquation.C

Contributors/Copyright:
    2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef forceEquation_H
#define forceEquation_H

#include "FieldValueExpressionDriver.H"

#include "fvMatrices/fvMatrix/fvMatrix.H"

#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceEquation Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class forceEquation
:
    protected FieldValueExpressionDriver
{
    // Private data

    forceEquation(const forceEquation&);

    exprString valueExpression_;

    exprString maskExpression_;

    bool verbose_;

    bool getMask(DynamicList<label> &,const word &psi);

public:

    // Constructors

        //- Construct from a dictionary
        forceEquation
        (
            const dictionary& ,
            const fvMesh&
        );

    // Destructor

        virtual ~forceEquation();

    //- fix equations
    void operator()(fvMatrix<T> &);

    //- where are the equations fixed
    tmp<volScalarField> getMask();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
