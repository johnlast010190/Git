/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2019-2025 Engys Ltd.

Class
    Foam::surfaceMap

Description
    Class designed to create a symmetric scalar field from an unsymmetric field.

SourceFiles
    surfaceMap.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMap_H
#define surfaceMap_H

#include "helyxMap.H"
#include "meshes/primitiveShapes/plane/plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class surfaceMap Declaration
\*---------------------------------------------------------------------------*/

class surfaceMap
:
    public helyxMap
{
protected:

    // Protected Data

        autoPtr<plane> mirror_;

        const fvMesh* mesh_;

        const Time* runTime_;

        bool imap_;
        bool bmap_;
        bool byType_;
        bool byPhysicalType_;
        bool byNameList_;

        word bodySurfaceType_;


public:

    // Public Data

        //- Patch names separated by comma
        std::vector<word> sourcePatchNames_;

        std::map<word, KNN> bndMaps_;

        //- Boundary tree ID => boundary region face ID maps
        std::map<word, std::vector<label>> faceIdMaps_;

        //- Boundary face centres
        std::map<word, DynamicList<point>> bSourceXYZ_;

        //- Boundary scalar field, bndName, scalarName: field
        std::map<word, DynamicList<scalar>> bSourceScalarFields_;


    // Constructors

        //- Default construct
        surfaceMap();

        //- Construct from fvMesh
        surfaceMap(const fvMesh& mesh);

        //- Construct from fvMesh and Time
        surfaceMap(const fvMesh* mesh, const Time* runTime);


    //- Destructor
    virtual ~surfaceMap();


    // Member Functions

        bool internalSym()
        {
            return imap_;
        }

        bool bndSym()
        {
            return bmap_;
        }

        void setInput(const dictionary& dict);

        void setInput(const dictionary& dict, const argList& args);

        void createMirrorFields(const word& mapTime);

        plane& mirror()
        {
            return mirror_();
        }

        label patchId(const word& name);

        word mirrorName(const word& name)
        {
            return name + "_symp";
        }

        void getParallelFields();

        void getMirrorFields();

        void getSymFields();

        void buildSearchTree();

        void buildSearchTreePar();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
