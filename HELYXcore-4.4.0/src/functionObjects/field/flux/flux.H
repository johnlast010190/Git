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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::functionObjects::flux

Group
    grpFieldFunctionObjects

Description
    Calculates the flux of a field.  The operation is limited to
    surfaceVectorFields and volVectorFields, and the output is a
    surfaceScalarField.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    flux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_flux_H
#define functionObjects_flux_H

#include "fieldExpression/fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class flux Declaration
\*---------------------------------------------------------------------------*/

class flux
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the flux of a volVectorField and register the result
        template<class FieldType, class RhoFieldType>
        bool calcVolFlux(const RhoFieldType& rho);

        //- Calculate the flux of a surfaceVectorField and register the result
        template<class FieldType, class RhoFieldType>
        bool calcSurFlux(const RhoFieldType& rho);

        //- Calculate the flux field and return true if successful
        virtual bool calc() override;

        //- rho flied name
        word rhoName_;


public:

    //- Runtime type information
    TypeName("flux");


    // Constructors

        //- Construct from Time and dictionary
        flux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flux();

    virtual bool execute() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flux/fluxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
