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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::thermoCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    thermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef thermoCloud_H
#define thermoCloud_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class thermoCloud Declaration
\*---------------------------------------------------------------------------*/

class thermoCloud
{
public:

    //- Runtime type information
    TypeName("thermoCloud");

    // Constructors

        //- Null constructor
        thermoCloud();

        //- Disallow default bitwise copy construction
        thermoCloud(const thermoCloud&) = delete;


    //- Destructor
    virtual ~thermoCloud();


    // Member Functions

        //- Return tmp equivalent particulate emission
        virtual tmp<volScalarField> Ep() const = 0;

        //- Return tmp equivalent particulate absorption
        virtual tmp<volScalarField> ap() const = 0;

        //- Return tmp equivalent particulate scattering factor
        virtual tmp<volScalarField> sigmap() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermoCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
