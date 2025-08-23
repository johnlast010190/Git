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

Class
    Foam::RepositoryBase

Description
    Base for the repositories. Manly here to generate a list so that
    all repositories can be updated in one go

SourceFiles
    RepositoryBase.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef RepositoryBase_H
#define RepositoryBase_H

#include "containers/Lists/DynamicList/DynamicList.H"
#include "db/objectRegistry/objectRegistry.H"

#include "helpers/DebugOStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RepositoryBase Declaration
\*---------------------------------------------------------------------------*/

class RepositoryBase
:
    public regIOobject
{
    // Private data

    // Constructors

    //- List of the repositories
    static DynamicList<RepositoryBase *> allRepos_;

protected:

    //- Construct null
    RepositoryBase(const IOobject &);

public:
     //- Runtime type information
    TypeName("RepositoryBase");

    //- Destructor

    virtual ~RepositoryBase();

    //- Update all the repositories
    static void updateRepos();

    //- method to update this repo for a new timestep
    virtual void updateRepo() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
