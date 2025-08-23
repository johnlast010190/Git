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
    (c) 2025 Engys Ltd.

Class
    Foam::functionObjects::runTimeVis::rtppFunctionObject

Group
    grpGraphicsFunctionObjects

Description
    Base class for all RTPP-based function objects

SourceFiles
    rtppFunctionObject.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "db/functionObjects/stateFunctionObject/stateFunctionObject.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam::functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class rtppFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class rtppFunctionObject : public stateFunctionObject
{
protected:
    // Constructors
    rtppFunctionObject
        (
            const word &name,
            const Time &runTime,
            const dictionary& dict
        );

    //- Read the post-processing controls
    bool read(const dictionary &dict) override;

    //- Destructor
    ~rtppFunctionObject() override;

protected:

    static bool readDebugFlag(const Foam::dictionary &dictionary);

    void finishInitialization() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
