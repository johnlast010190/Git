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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "simpleDataFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"
#include "include/OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleDataFunctionObject, 0);

fileName simpleDataFunctionObject::defaultPostProcDir_("postProcessing");

void simpleDataFunctionObject::setPostProcDir(const fileName &f)
{
    Info<< "Setting output directory name for simpleFunctionObjects to "
        << f << endl;
    defaultPostProcDir_=f;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

simpleDataFunctionObject::simpleDataFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict),
    postProcDir_(defaultPostProcDir_)
{
    Dbug<< name << " - Constructor" << endl;

    if (dict.found("postProcDir")) {
        postProcDir_=fileName(
            dict.lookup("postProcDir")
        );
        Info<< name << " writes to " << postProcDir_
            << " instead of " << defaultPostProcDir_ << endl;
    }
}

fileName simpleDataFunctionObject::dataDir()
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    // make sure that when starting we take the start time
    if (
        obr_.time().timeIndex()
        <=
        obr_.time().startTimeIndex()+1
    ) {
        return baseDir()/obr_.time().timeName(
            obr_.time().startTime().value()
        );
    } else {
        return baseDir()/obr_.time().timeName();
    }
#else
    return baseDir()/obr_.time().timeName();
#endif
}

fileName simpleDataFunctionObject::baseDir()
{
    fileName theDir;
    fileName dir=dirName()+"_"+name();
    if (obr().name()!=polyMesh::defaultRegion) {
        dir=obr().name() / dir;
    }
    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        theDir =
            obr_.time().path()
            /".."
            /postProcDir_
            /dir;
    }
    else
    {
        theDir =
            obr_.time().path()
            /postProcDir_
            /dir;
    }

    return theDir;
}

bool simpleDataFunctionObject::start()
{
    Dbug<< name() << "::start()" << endl;

    simpleFunctionObject::start();

    mkDir(dataDir());

    return true;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
