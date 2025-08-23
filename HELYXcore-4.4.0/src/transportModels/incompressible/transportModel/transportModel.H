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
    (c) 2010-2024 Engys Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::transportModel

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    transportModel.C

\*---------------------------------------------------------------------------*/

#ifndef transportModel_H
#define transportModel_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "fields/volFields/volFieldsFwd.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class transportModel Declaration
\*---------------------------------------------------------------------------*/

class transportModel
{
public:

    //- Runtime type information
    TypeName("transportModel");


    // Constructors

        //- Null construct
        transportModel();

        //- Disallow default bitwise copy construction
        transportModel(const transportModel&) = delete;


    //- Destructor
    virtual ~transportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;

        //- Return constant density
        virtual tmp<volScalarField> rho() const = 0;

        //- Constant heat capacity
        virtual const volScalarField& Cp() const = 0;

        //- Constant thermal conductivity
        virtual tmp<volScalarField> lambda() const = 0;

        //- Constant turbulent Prandtl number
        virtual tmp<volScalarField> Prt() const = 0;

        //- Laminar Prandtl number
        virtual tmp<volScalarField> Pr() const = 0;

        //- Constant turbulent Prandtl number
        virtual tmp<volScalarField> alphaLam() const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
