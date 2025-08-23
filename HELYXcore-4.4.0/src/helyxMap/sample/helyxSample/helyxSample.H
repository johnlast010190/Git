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
    Foam::helyxSample

Description
    Class designed to sample an arbitrary number of flow fields from
    an unstructured mesh region into a uniform grid,
    or to another unstructured mesh region.

SourceFiles
    helyxSample.C

\*---------------------------------------------------------------------------*/

#ifndef helyxSample_H
#define helyxSample_H

#include <sys/stat.h>
#include "helyxMap.H"
#include "sample/gridField/gridField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class helyxSample Declaration
\*---------------------------------------------------------------------------*/

class helyxSample
:
    public helyxMap
{
    // Private Member Functions

        //- Helper function to initialise defaults
        void init();


protected:

    // Protected Data

        autoPtr<boundBox> sampleBox_;

        word databasePath_;


public:

    // Static Data Members

        static word databaseDir_;


    // Public Data

        //- Boundary map for the bodies
        KNN bodyBoundaryTree_;

        const fvMesh* mesh_;

        const Time* runTime_;

        std::map<word, gridField> sampleGrids_;

        bool toDatabase_;

        word rootCase_;

        bool visualize_;
        bool writeSolidMask_;
        bool writeSamples_;

        label nProc_;

        DynamicList<std::pair<label, label>> bodyFaces_;
        DynamicList<std::pair<label, label>> procIfaces_;

        word bodySurfaceType_;

        bool byPhysicalType_;
        bool byType_;
        bool byName_;
        bool byNameList_;

        label nx_;
        label ny_;
        label nz_;

        wordList bodyPatchNames_;


    // Static Member Functions

        //- Write HELYX banner to output stream
        static void writeBanner
        (
            Ostream& os,
            const word classType,
            const word objectName
        );


    // Constructors

        //- Default construct
        helyxSample();

        //- Construct from Time and fvMesh
        helyxSample(const Time* runTime, const fvMesh* mesh);


    //- Destructor
    virtual ~helyxSample();


    // Member Functions

        boundBox& sampleBox()
        {
            return sampleBox_();
        }

        virtual void setOptions(const argList& args);

        bool parRun()
        {
            return (nProc_ > 1);
        }

        bool isProcInterface(const word& name)
        {
            word wd = "procBoundary";
            label j = name.find(wd);

            return (j >= 0);
        }

        virtual void getMapTime()
        {}

        virtual void combineFields()
        {}

        virtual void constructKnn()
        {}

        virtual void getBodySurface()
        {}

        void getBodyPatchNames();

        virtual void setInput(const dictionary& dict, const argList& args);

        word databasePath() const
        {
            return databasePath_;
        }

        void getSampleData(const word& funct = "field");

        void buildSearchTrees();

        virtual void searchNodeField
        (
            gridField& grid,
            const word& funct = "field"
        );

        void getGeomSample(gridField& grid);

        void writeNodeField(gridField& grid);

        void saveObjectMask(gridField& grid, const word& type);

        void saveObjectMasks(const word& type);

        //- Save block samples to database
        void saveToDatabase(const word& funct = "field");

        void saveToDatabase(gridField& grid, const word& funct = "field");

        void writeNodeField();

        void writeVisualization();

        void writeMeshDict(gridField& grid);

        void writeScalarField(gridField& grid, label k = 0);

        void writeSolidMask(gridField& grid);

        void writeVectorField(gridField& grid, label k = 0);

        //- Get the scalar field at a nodal point by kdTree search
        inline void getGridField
        (
            gridField& grid,
            label i,
            label j,
            label k,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist,
            bool insolid
        );

        inline void outDomainGrid(gridField& grid, label i, label j, label k);

        //- Get the vector field at a nodal point
        inline void getGridVectorField
        (
            gridField& grid,
            label i,
            label j,
            label k,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist,
            bool insolid
        );

        inline scalar interp
        (
            const word& name,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );

        inline point vInterp
        (
            const word& name,
            const std::vector<label>& nb,
            const std::vector<scalar>& dist
        );
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline bool isFile(const char* path)
{
    struct stat buf;
    stat(path, &buf);

    return S_ISREG(buf.st_mode);
}

inline bool isDir(const char* path)
{
    struct stat buf;
    stat(path, &buf);

    return S_ISDIR(buf.st_mode);
}

//- Get files or sub-directories
void getFiles
(
    wordList& v,
    const word& root = "./",
    const word& contains = "none"
);

void readDirectory
(
    wordList& v,
    const word& root = "./",
    const word& filter = "none"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helyxSampleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
