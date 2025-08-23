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
    Foam::helyxMap

Description
    A set of classes designed to map an arbitrary number of flow fields from
    an unstructured mesh region into a uniform grid, or to another unstructured
    mesh region.

SourceFiles
    helyxMap.C

\*---------------------------------------------------------------------------*/

#ifndef helyxMap_H
#define helyxMap_H

#include "tools/readFields/readFields.H"
#include "memInfo/memInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class helyxMap Declaration
\*---------------------------------------------------------------------------*/

class helyxMap
{
protected:

    // Protected Data

        bool norm_;

        label parModel_;

        scalar extCoef_;

        labelHashSet wallDistPatchs_;

        wordReList wallPatchNames_;

        word bboxType_;

        label boxCellNum_;

        //- Report memory usage for debug
        bool reportMem_;

        bool scaleSource_;


public:

    // Public Data

        typedef Tuple2<label, label> index;

        word caseDir_;
        word caseType_;

        label nProcs_;
        label myId_;
        label masterId_;

        scalar deltaBox_;
        DynamicList<label> procs_;

        bool wDistMap_;

        DynamicList<bool> cellInDomain_;

        DynamicList<point> sourceXyz_;
        std::map<word, DynamicList<scalar>> sourceScalarFields_;
        std::map<word, DynamicList<point>> sourceVectorFields_;
        DynamicList<scalar> sourceWDists_;
        std::map<label, label> cellMap_;
        std::vector<KNN> kdTrees_;
        KNN internalMap_;
        scalar deltaY_;

        scalar sourceAlphaMax_;

        //- Parallel run alphaMap
        Array2d<label> alphaMap_;

        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;

        DynamicList<label> toProcs_;
        DynamicList<label> fromProcs_;

        DynamicList<index> interProcs_;
        scalar alphaMax_;
        scalar errorBound_;

        scalar rhoRefSource_;
        scalar rhoRefTarget_;

        //- Rotating degree from source to target, anti-clockwise positive
        scalar UrotDegreeFromSource_;

        vector UrefSource_;
        vector UrefTarget_;
        bool interp_;
        bool mapBoundary_;
        wordList mapFixedBCs_;
        word function_;
        word mapTimeName_;
        word tgtTimeName_;
        HashTable<word> fieldTypes_;

        word sourceCase_;
        HashSet<word> mapScalarFields_;
        HashSet<word> mapVectorFields_;
        wordList mapFieldNames_;

        HashSet<word> mapSurfaceScalarFields_;
        label nwdist_;
        word mapCase_;

        wordList sourceRegions_;
        wordList targetRegions_;
        wordList targetRegionDirs_;
        bool allRegions_;

        autoPtr<class readFields> source_;
        autoPtr<class readFields> target_;


    // Static Member Functions

        static label trimLeft(word& str, const word& del = " ");

        static label trimRight(word& str, const word& del = " ");

        static void trim(word& str, const word& del = " ");

        static void replaceString(word& str);

        static word intToString(label i)
        {
            std::stringstream ss;
            ss << i;
            return ss.str();
        }

        static word doubleToString(scalar d)
        {
            std::stringstream ss;
            ss << d;
            return ss.str();
        }

        static void stringTokenize
        (
            wordList& container,
            const word& in,
            const word& delimiters
        );


    // Constructors

        //- Default construct
        helyxMap();


    //- Destructor
    virtual ~helyxMap();


    // Member Functions

        // Access

            const labelHashSet& wallDistPatchs()
            {
                return wallDistPatchs_;
            }

            bool reportMemUsage()
            {
                return reportMem_;
            }

            bool scaleSource()
            {
                return scaleSource_;
            }

            label parModel()
            {
                return parModel_;
            }

            scalar extCoef()
            {
                return extCoef_;
            }

        void reportMemory(memInfo& mem, const word& desc, scalar& maxUsage);

        void getWallDistPatchs(const fvMesh& mesh);

        bool inList(label k, const DynamicList<label>& list);

        bool inList(const word& key, const wordList& list);

        //- Build parallel run kdTrees
        virtual void buildKdTreesPar();

        virtual void constructInternalKnn(scalar cutoff = 0.8);

        scalar xbar(scalar x)
        {
            return target().xbar(x);
        }

        scalar ybar(scalar y)
        {
            return target().ybar(y);
        }

        scalar zbar(scalar z)
        {
            return target().zbar(z);
        }

        label inGroup(label i)
        {
            if (i > sourceWDists_.size() - 1)
            {
                i = sourceWDists_.size() - 1;
            }

            return label(sourceWDists_[i]/deltaY_);
        }

        bool master()
        {
            return (myId_ == masterId_);
        }

        label wait()
        {
            if (!parRun())
            {
                return 0;
            }
            else
            {
                return gMax(procs_);
            }
        }

        word caseName()
        {
            if (!parRun())
            {
                wordList vtmp;
                stringTokenize(vtmp, caseDir_, "/");

                word cname = "fields";
                if (vtmp.size() == 1)
                {
                    cname += caseDir_;
                }
                else
                {
                    for (label i = 0; i < vtmp.size(); i++)
                    {
                        cname += ("_" + vtmp[i]);
                    }
                }
                return cname;
            }

            word ss = "processor";
            ss += intToString(myId_);
            return ss;
        }

        void parInit(label np);

        //- For parallel runs, get domain fields from source case
        void getDomainFields();

        //- Single-field parallel model
        void getDomainFields0();

        void createSource(const fvMesh* mesh, const Time* runTime);

        void createTarget(const fvMesh* mesh, const Time* runTime);

        void clearSourceTarget();

        label searchError();

        //- Mapping error based on a scalar
        scalar mapError(const word& scname);

        scalar UmapError(const word& vecname);

        void mapScalarField(const word& scname);
        void mapVectorField(const word& vecname);
        void mapScalarFieldPar(const word& scname);
        void mapVectorFieldPar(const word& vecname);

        void getFieldTypes(const fvMesh& mesh, const word& timeName);

        void scalePressure
        (
            const dimensionSet& sourceDim,
            const dimensionSet& targetDim
        );

        void checkFieldDimensions();

        void mapFields();

        void setMapTime(Time& runTime);

        virtual void setInput(const dictionary& dict);

        virtual void setInput(const dictionary& dict, const argList& args);

        bool parRun()
        {
            return (nProcs_ > 1);
        }

        class readFields& source()
        {
            return source_();
        }

        class readFields& target()
        {
            return target_();
        }

        virtual void setOptions(const argList& args);

        label cellIndex
        (
            scalar alphai,
            label kdTreeId
        )
        {
            if (alphai > sourceAlphaMax_ || !wDistMap_)
            {
                return cellMap_[kdTreeId];
            }
            else
            {
                label ialpha = int(alphai/deltaY_);
                return alphaMap_(ialpha, kdTreeId);
            }
        }

        bool inDomain(const point& pt)
        {
            return target().inDomain(pt);
        }

        inline void searchNearest
        (
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar x,
            scalar y,
            scalar z
        );

        inline void searchNearest
        (
            std::vector<label>& nb,
            std::vector<scalar>& dist,
            scalar alphai,
            scalar x,
            scalar y,
            scalar z
        );
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
void interProcTrans
(
    DynamicList<T>& field,
    const DynamicList<point>& verts,
    const List<boundBox>& tboxs,
    const List<boundBox>& sboxs,
    const DynamicList<T>& source,
    const word& name
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helyxMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
