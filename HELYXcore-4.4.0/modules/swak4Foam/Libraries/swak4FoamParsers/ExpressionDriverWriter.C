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

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "ExpressionDriverWriter.H"

#include "CommonValueExpressionDriver.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ExpressionDriverWriter, 0);

ExpressionDriverWriter::ExpressionDriverWriter(
        const word &name,
        CommonValueExpressionDriver &driver
):
    regIOobject(
        IOobject(
            name,
            driver.mesh().time().timeName(),
            "swak4Foam",
            driver.mesh().time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    driver_(driver)
{
    if (debug) {
        Pout<< "ExpressionDriverWriter at " << objectPath() << " created" << endl;
    }

    if (headerOk()) {
        if (debug) {
            Pout<< "Found a file " <<  objectPath() << endl;
        }

        readData(
            readStream(
                "ExpressionDriverWriter"
#ifdef FOAM_READSTREAM_METHOD_NEEDS_BOOL_PARAMETER
                ,true
#endif
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ExpressionDriverWriter::~ExpressionDriverWriter()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ExpressionDriverWriter::writeData(Ostream &os) const
{
    if (debug) {
        Pout<< "ExpressionDriverWriter at " << objectPath()
            << " writing" << endl;
    }

    dictionary dict;

    driver_.prepareData(dict);

    dict.write(os,false);

    if (debug) {
        Pout<< "written " << dict << endl;
    }

    return os.good();
}

bool ExpressionDriverWriter::readData(Istream &is)
{
    if (debug) {
        Pout<< "ExpressionDriverWriter at " << objectPath()
            << " reading" << endl;
    }

    const dictionary dict(is);

    if (debug) {
        Pout<< "reading " << dict << endl;
    }

    driver_.getData(dict);

    return !is.bad();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
