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
    (c) 2024 Engys Ltd

Description
    Macro definitions for debug switches.

\*---------------------------------------------------------------------------*/

#ifndef defineDebugSwitch_H
#define defineDebugSwitch_H

#include "global/debug/simpleRegIOobject.H"
#include "global/debug/debug.H"
#include "primitives/ints/label/label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class RegisterDebugSwitchImpl final
:
    public ::Foam::simpleRegIOobject
{
    int* flagAddr;
public:

    RegisterDebugSwitchImpl(int* flagAddress, const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addDebugObject, name),
        flagAddr(flagAddress)
    {}

    virtual ~RegisterDebugSwitchImpl() = default;

    virtual void readData(Foam::Istream& is)
    {
        *flagAddr = readLabel(is);
    }

    virtual void writeData(Foam::Ostream& os) const
    {
        os << *flagAddr;
    }
};

//- Define the debug information, lookup as \a Name
template<class Type>
class RegisterDebugSwitch final {
public:
    //- The unique RegisterDebugSwitch object
    static const RegisterDebugSwitchImpl registerDebugSwitch;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Define the debug information, lookup as \a Name
#define registerTemplateDebugSwitchWithName(Type,Name)                         \
    template<>                                                                 \
    const Foam::RegisterDebugSwitchImpl                                        \
        Foam::RegisterDebugSwitch<Type>::registerDebugSwitch(&Type::debug, Name); \


//- Define the debug information, lookup as \a Name
#define registerDebugSwitchWithName(Type,Tag,Name)                             \
    static RegisterDebugSwitchImpl add##Tag##ToDebug_(&Type::debug, Name)


//- Define the debug information, lookup as \a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch)                     \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch))

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch)                                   \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch);           \
    registerDebugSwitchWithName(Type, Type, Type::typeName_())

//- Define the debug information for templates, lookup as \a Name
#define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)             \
    template<>                                                                 \
    defineDebugSwitchWithName(Type, Name, DebugSwitch);                        \
    registerTemplateDebugSwitchWithName(Type, Name)

//- Define the debug information for templates sub-classes, lookup as \a Name
#define defineTemplate2DebugSwitchWithName(Type, Name, DebugSwitch)            \
    template<>                                                                 \
    defineDebugSwitchWithName(Type, Name, DebugSwitch);                        \
    registerTemplateDebugSwitchWithName(Type, Name)

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, DebugSwitch)                           \
    defineTemplateDebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, DebugSwitch)                      \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
