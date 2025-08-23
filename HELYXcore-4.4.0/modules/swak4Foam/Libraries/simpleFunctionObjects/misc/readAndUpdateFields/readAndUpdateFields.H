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

Class
    Foam::readAndUpdateFields

Description
    Reads fields at startup and updates the boundary conditions at every
    timestep. Also makes sure that the fields are automatically written

SourceFiles
    readAndUpdateFields.C

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef readAndUpdateFields_H
#define readAndUpdateFields_H

#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "include/swak.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyTopoChangeMap;
class polyMeshMap;
class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class readAndUpdateFields Declaration
\*---------------------------------------------------------------------------*/

class readAndUpdateFields
{
protected:

    // Protected data

        //- Name of this set of readAndUpdateFields object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Fields to load
        wordList fieldSet_;

        //- Loaded fields
        PtrList<volScalarField> vsf_;
        PtrList<volVectorField> vvf_;
        PtrList<volSphericalTensorField> vSpheretf_;
        PtrList<volSymmTensorField> vSymmtf_;
        PtrList<volTensorField> vtf_;

        PtrList<pointScalarField> psf_;
        PtrList<pointVectorField> pvf_;
        PtrList<pointSphericalTensorField> pSpheretf_;
        PtrList<pointSymmTensorField> pSymmtf_;
        PtrList<pointTensorField> ptf_;

        PtrList<surfaceScalarField> ssf_;
        PtrList<surfaceVectorField> svf_;
        PtrList<surfaceSphericalTensorField> sSpheretf_;
        PtrList<surfaceSymmTensorField> sSymmtf_;
        PtrList<surfaceTensorField> stf_;

        autoPtr<pointMesh> pMesh_;

        //- Should the boundary conditions be corrected at every time-step?
        bool correctBoundary_;


    // Protected Member Functions

        const pointMesh& pMesh(const polyMesh& mesh);

        template<class Type>
        bool loadField
        (
            const word&,
            PtrList<VolField<Type>>&,
            PtrList<SurfaceField<Type>>&,
            PtrList<PointField<Type>>&
        );

        template <class FType>
        void correctBoundaryConditions(PtrList<FType>&);


public:

    //- Runtime type information
    TypeName("readAndUpdateFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        readAndUpdateFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        readAndUpdateFields(const readAndUpdateFields&) = delete;


    //- Destructor
    virtual ~readAndUpdateFields();


    // Member Functions

        //- Return name of the readAndUpdateFields object
        virtual const word& name() const
        {
            return name_;
        }

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();
        virtual void setPostProcess(const bool&){};

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
        virtual bool write();
#else
        virtual void write();
#endif

        //- Update for changes of mesh
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update for changes of mesh
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Update for changes of mesh
#ifdef FOAM_MOVEPOINTS_GETS_POLYMESH
        virtual void movePoints(const polyMesh&)
#else
        virtual void movePoints(const pointField&)
#endif
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const readAndUpdateFields&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "readAndUpdateFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
