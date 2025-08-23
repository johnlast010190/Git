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
    Foam::SurfacesRepository

Description

SourceFiles
    SurfacesRepositoryI.H
    SurfacesRepository.C
    SurfacesRepositoryIO.C

Contributors/Copyright:
    2012-2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SurfacesRepository_H
#define SurfacesRepository_H

#include "containers/HashTables/PtrMap/PtrMap.H"
#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "RepositoryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SurfacesRepository Declaration
\*---------------------------------------------------------------------------*/

class SurfacesRepository
:
    public RepositoryBase
{
    // Private data

    //- the actual surfaces
    HashPtrTable<sampledSurface,word> surfaces_;

    //- names of the formats to be written
    HashTable<word,word> formatNames_;

    // Constructors

    //- Construct null
    SurfacesRepository(const IOobject &);

    //- The only instance of this repository
    static SurfacesRepository *repositoryInstance;

public:
     //- Runtime type information
    TypeName("SurfacesRepository");

    //- Destructor

    virtual ~SurfacesRepository();

    //- Get the Singleton of this class
    static SurfacesRepository &getRepository(const objectRegistry &mesh);

    //- get a surface using a dictionary (surface will be constructed if not existing)
    sampledSurface &getSurface(
        const dictionary &dict,
        const fvMesh &mesh
    );

    //- get a surface using its name (surface must exist)
    sampledSurface &getSurface(
        const word &name,
        const fvMesh &mesh
    );

    //- force an update of a surface
    bool updateSurface(
        const word &name,
        const fvMesh &mesh
    );

    //- Make sure that this class is not abstract
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
