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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "RepositoryBase.H"

#include "include/swakTime.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

defineTypeNameAndDebug(RepositoryBase, 0);

DynamicList<RepositoryBase *> RepositoryBase::allRepos_;

RepositoryBase::RepositoryBase(const IOobject &o)
    :
    regIOobject(o)
{
    Dbug<< "Adding Repo " << getHex(this) << " to global chain" << endl;

    allRepos_.append(this);
    allRepos_.shrink();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

RepositoryBase::~RepositoryBase()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void RepositoryBase::updateRepos()
{
    Sbug << "Updating all repositories" << endl;

    typedef DynamicList<RepositoryBase *> repoList;

    forAllIter(repoList,allRepos_,it)
    {
    (*it)->updateRepo();
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace

// ************************************************************************* //
