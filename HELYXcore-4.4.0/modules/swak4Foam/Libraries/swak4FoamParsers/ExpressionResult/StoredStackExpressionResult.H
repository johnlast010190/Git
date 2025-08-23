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
    (c) 1991-2008 OpenCFD Ltd.

Class
    Foam::StoredStackExpressionResult

Description

SourceFiles
    StoredStackExpressionResultI.H
    StoredStackExpressionResult.C
    StoredStackExpressionResultIO.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef StoredStackExpressionResult_H
#define StoredStackExpressionResult_H

#include "StackExpressionResult.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class StoredStackExpressionResult Declaration
\*---------------------------------------------------------------------------*/

class StoredStackExpressionResult
:
    public StackExpressionResult
{
    // Private data

protected:

    //- reset at new timestep according to type
    virtual void resetInternal();

public:

     //- Runtime type information
    TypeName("StoredStackExpressionResult");

    // Constructors

    //- Construct null
    StoredStackExpressionResult();

    //- Construct as copy
    StoredStackExpressionResult(const StoredStackExpressionResult&);

    //- Construct from a dictionary
    StoredStackExpressionResult(const dictionary&);

    // Destructor

    virtual ~StoredStackExpressionResult();

    void operator=(const StoredStackExpressionResult&);

    virtual autoPtr<ExpressionResult> clone() const {
        return autoPtr<ExpressionResult>(
            new StoredStackExpressionResult(*this)
        );
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StoredStackExpressionResultI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
