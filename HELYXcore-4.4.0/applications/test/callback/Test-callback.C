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
    (c) 2011-2015 OpenFOAM Foundation

Application
    callBackTest

Description

\*---------------------------------------------------------------------------*/

#include "Callback.H"

using namespace Foam;

class callback
:
    public Callback<callback>
{
public:

    callback(CallbackRegistry<callback>& cbr)
    :
        Callback<callback>(cbr)
    {}

    ~callback()
    {}

    virtual const word& name() const
    {
        return word::null;
    }

    void testCallbackFunction() const
    {
        Info<< "calling testCallbackFunction for object " << name() << endl;
    }
};


class callbackRegistry
:
    public CallbackRegistry<callback>
{
public:

    callbackRegistry()
    {}

    ~callbackRegistry()
    {}

    void testCallbackFunction() const
    {
        forAllConstIter(callbackRegistry, *this, iter)
        {
            iter().testCallbackFunction();
        }
    }
};


class objectWithCallback
:
    public callback
{
    word name_;

public:

    objectWithCallback(const word& n, callbackRegistry& cbr)
    :
        callback(cbr),
        name_(n)
    {}

    virtual const word& name() const
    {
        return name_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    callbackRegistry cbr;

    objectWithCallback ob1("ob1", cbr);
    objectWithCallback ob2("ob2", cbr);

    cbr.testCallbackFunction();

    {
        objectWithCallback ob1("ob1", cbr);
        cbr.testCallbackFunction();
    }

    cbr.testCallbackFunction();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
