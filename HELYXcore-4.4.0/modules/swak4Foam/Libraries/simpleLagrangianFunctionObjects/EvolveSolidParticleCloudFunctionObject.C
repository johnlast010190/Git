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
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "EvolveSolidParticleCloudFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    // Specialization because solidParticleCloud has no evolve
template<>
bool EvolveCloudFunctionObject<solidParticleCloud>::execute(bool forceWrite)
{
#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    if (!cloud_.valid()) {
        this->start();
    }

    if (
        lastTimeStepExecute_
        !=
        obr().time().timeIndex()
    ) {
        lastTimeStepExecute_=obr().time().timeIndex();
    } else {
        return false;
    }
#endif

    Info<< "Moving solidParticeCloud:" << cloud_->name()
        << " with " << cloud_->size() << " particles" << endl;
    cloud_->move(g());
    Info<< tab << cloud_->size() << " particles after moving"
        << endl;

    if (
        obr().time().outputTime()
        ||
        forceWrite
    ) {
        Info<< "Writing cloud " << cloud_->name() << endl;
        cloud_->write();
    }

    return true;
}


    defineTypeNameAndDebug(EvolveSolidParticleCloudFunctionObject, 0);

    addNamedToRunTimeSelectionTable
    (
        functionObject,
        EvolveSolidParticleCloudFunctionObject,
        dictionary,
        evolveSolidParticleCloud
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

EvolveSolidParticleCloudFunctionObject::EvolveSolidParticleCloudFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    EvolveCloudFunctionObject<solidParticleCloud>(
        name,
        t,
        dict
    )
{
    Dbug<< this->name() << " Construktor" << endl;

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
    this->read(dict);
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool EvolveSolidParticleCloudFunctionObject::start()
{
    Dbug<< this->name() << "::start()" << endl;

    cloud().set(
        new solidParticleCloud(
            //            dynamicCast<const fvMesh &>(
            dynamic_cast<const fvMesh &>(
                obr()
            ),
            cloudName()
        )
    );

    return true;
}


} // namespace Foam

// ************************************************************************* //
