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

#include "updateSimpleFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(updateSimpleFunctionObject, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

updateSimpleFunctionObject::updateSimpleFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    simpleFunctionObject(name,t,dict)
{
    Pout<< "Creating " << name << endl;

    runIfStartTime_=dict_.lookupOrDefault<bool>("runIfStartTime",false);
    onlyAtStartup_=dict_.lookup<bool>("onlyAtStartup");
    if (onlyAtStartup_) {
        runIfStartTime_=dict_.lookup<bool>("runIfStartTime");
    }
}

bool updateSimpleFunctionObject::start()
{
    Pbug << "start() started" << endl;
    simpleFunctionObject::start();
    Pbug << "start() called parent" << endl;

    if (onlyAtStartup_) {
        if (
            !runIfStartTime_
            ||
            time().timeIndex()==0
        ) {
            Pbug << "Calling recalc()" << endl;
            recalc();
        }
    }

    Pbug << "start() ended" << endl;

    return true;
}

void updateSimpleFunctionObject::writeSimple()
{
    Pbug << "write() started" << endl;

    if (!onlyAtStartup_) {
        Pbug << "Calling recalc() always" << endl;
        recalc();
    }

    Pbug << "write() ended" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
