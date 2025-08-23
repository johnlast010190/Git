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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::functionObjects::scale

Group
    grpFieldFunctionObjects

Description
    Calculates the re-scaling of a field.  The operation is limited to
    surfaceScalarFields and volVectorFields, and the output is a volScalarField.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scale.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scale_H
#define functionObjects_scale_H

#include "fieldExpression/fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class scale Declaration
\*---------------------------------------------------------------------------*/

class scale
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the re-scaling of either a
        //  volScalarField or a surfaceScalarField and register the result
        template<class FieldType>
        bool calcScale();

        //- Calculate the scaled field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("scale");


    // Constructors

        //- Construct from Time and dictionary
        scale
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~scale();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "scale/scaleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
