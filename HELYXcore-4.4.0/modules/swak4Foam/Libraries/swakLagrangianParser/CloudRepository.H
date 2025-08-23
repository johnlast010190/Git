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
    Foam::CloudsRepository

Description
    Store clouds that nobody else wants

SourceFiles
    CloudRepository.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CloudRepository_H
#define CloudRepository_H

#include "containers/HashTables/PtrMap/PtrMap.H"
#include "fields/cloud/cloud.H"
#include "repositories/RepositoryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class ReaderParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class CloudRepository Declaration
\*---------------------------------------------------------------------------*/

class CloudRepository
:
    public RepositoryBase
{
    // Private data

    //- the actual sets
    HashPtrTable<cloud,word> clouds_;

    //- Clouds that can be updated
    HashPtrTable<ReaderParticleCloud,word> updateableClouds_;

   // Constructors

    //- Construct null
    CloudRepository(const IOobject &);

    //- The only instance of this repository
    static CloudRepository *repositoryInstance;

public:
     //- Runtime type information
    TypeName("CloudRepository");

    //- Destructor

    virtual ~CloudRepository();

    //- Get the Singleton of this class
    static CloudRepository &getRepository(const objectRegistry &mesh);

    void addCloud(
        autoPtr<cloud> c
    );

    void addUpdateableCloud(
        autoPtr<ReaderParticleCloud> c
    );

    //- Make sure that this class is not abstract. Write sets
    bool writeData(Ostream &) const;

    //- update repo at the start of a time-step
    void updateRepo();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
