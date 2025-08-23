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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CloudRepository.H"

#include "meshes/polyMesh/polyMesh.H"
#include "meshSearch/meshSearch.H"

#include "ReaderParticle/ReaderParticleCloud.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

defineTypeNameAndDebug(CloudRepository, 0);

CloudRepository *CloudRepository::repositoryInstance(NULL);

CloudRepository::CloudRepository(const IOobject &o)
    :
    RepositoryBase(o)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CloudRepository::~CloudRepository()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

CloudRepository &CloudRepository::getRepository(const objectRegistry &obr)
{
    CloudRepository*  ptr=repositoryInstance;

    if (debug) {
        Pout<< "CloudRepository: asking for Singleton" << endl;
    }

    if (ptr==NULL) {
        Pout<< "swak4Foam: Allocating new repository for sampledSets\n";

        ptr=new CloudRepository(
            IOobject(
                "swakClouds",
                obr.time().timeName(),
                "cloudRepository",
                obr.time(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
    }

    repositoryInstance=ptr;

    return *repositoryInstance;
}

bool CloudRepository::writeData(Ostream &f) const
{
    if (debug) {
        Info<< "CloudRepository::write()" << endl;
    }

    f << "Non-updateable " << clouds_.toc() << endl;

    f << "Updateable " << updateableClouds_.toc() << endl;

    return true;
}

void CloudRepository::addUpdateableCloud(
    autoPtr<ReaderParticleCloud> c
) {
    const word &name=c->name();

    if (clouds_.found(name)) {
        FatalErrorIn("CloudRepository::addCloud")
            << "There is already a cloud " << name
                << " in the non-updateable clouds. I guess there is a mistake"
                << endl
                << exit(FatalError);
    }

    if (updateableClouds_.found(name)) {
        FatalErrorIn("CloudRepository::addCloud")
            << "Repository of updateable clouds already has an entry "
                << name << ". This can't be right"
                << endl
                << exit(FatalError);
    } else {
        updateableClouds_.insert(
            name,
            c.ptr()
        );
    }
}

void CloudRepository::addCloud(
    autoPtr<cloud> c
) {
    const word &name=c->name();

    if (updateableClouds_.found(name)) {
        FatalErrorIn("CloudRepository::addCloud")
            << "There is already a cloud " << name
                << " in the updateable clouds. I guess there is a mistake"
                << endl
                << exit(FatalError);
    }

    if (clouds_.found(name)) {
        WarningIn("CloudRepository::addCloud")
            << "Repository of clouds already has an entry "
                << name <<". Overwriting. Expect strange behaviour"
                << endl;
        clouds_.set(
            name,
            c.ptr()
        );
    } else {
        clouds_.insert(
            name,
            c.ptr()
        );
    }
}

void CloudRepository::updateRepo()
{
    clouds_.clear();

    typedef HashPtrTable<ReaderParticleCloud,word> updateTable;

    forAllIter(updateTable,updateableClouds_,it)
    {
        //    (*it)->clear();
        //        Info<< "Updating " << it.key() << endl;
        (*it)->clearData();
        ReaderParticle::readFields(*(*it));
    }

}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
