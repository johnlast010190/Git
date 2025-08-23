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
    (c) 2024 Engys Ltd.

Class
    Foam::timeSinceChangePluginFunction

Description
    Time since the state machine last changed its state

SourceFiles
    timeSinceChangePluginFunction.C

Contributors/Copyright:
    2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef timeSinceChangePluginFunction_H
#define timeSinceChangePluginFunction_H

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class timeSinceChangePluginFunction Declaration
\*---------------------------------------------------------------------------*/

template <class PluginType>
class timeSinceChangePluginFunction
:
    public PluginType
{
    word machineName_;


protected:

    void setArgument(label index, const word& value);


public:

    //- Runtime type information
    TypeName("timeSinceChangePluginFunction");


    // Constructors

        //- Construct from components
        timeSinceChangePluginFunction
        (
            const typename PluginType::driverType &parentDriver,
            const word &name
        );

        //- Disallow default bitwise copy construction
        timeSinceChangePluginFunction
        (
            const timeSinceChangePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~timeSinceChangePluginFunction() {}


    // Member Functions

        void doEvaluation();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const timeSinceChangePluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
