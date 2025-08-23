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
    Foam::isStatePluginFunction

Description
    Check if a state machine is in a certain state and return true if it is

SourceFiles
    isStatePluginFunction.C

Contributors/Copyright:
    2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef isStatePluginFunction_H
#define isStatePluginFunction_H

#include "plugins/FieldValuePluginFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class isStatePluginFunction Declaration
\*---------------------------------------------------------------------------*/

template <class PluginType>
class isStatePluginFunction
:
    public PluginType
{

    word machineName_;

    word state_;


protected:

    void setArgument(label index, const word& value);


public:

    //- Runtime type information
    TypeName("isStatePluginFunction");


    // Constructors

        //- Construct from components
        isStatePluginFunction
        (
            const typename PluginType::driverType &parentDriver,
            const word &name
        );

        //- Disallow default bitwise assignment
        isStatePluginFunction(const isStatePluginFunction&) = delete;


    //- Destructor
    virtual ~isStatePluginFunction() {}


    // Member Functions

        void doEvaluation();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const isStatePluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
