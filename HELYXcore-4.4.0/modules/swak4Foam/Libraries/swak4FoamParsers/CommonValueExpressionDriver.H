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

Class
    Foam::CommonValueExpressionDriver

Description

SourceFiles
    CommonValueExpressionDriverI.H
    CommonValueExpressionDriver.C

Contributors/Copyright:
    2010-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef CommonValueExpressionDriver_H
#define CommonValueExpressionDriver_H

#include "include/swak.H"

#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "interpolations/interpolationTable/interpolationTable.H"
#ifdef FOAM_HAS_INTERPOLATION2DTABLE
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
#endif
#ifndef FOAM_DLLIBRARY_USES_STATIC_METHODS
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#endif

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helpers/exprString.H"
#include "ExpressionResult/ExpressionResult.H"
#include "ExpressionResult/StoredExpressionResult.H"
#include "ExpressionResult/DelayedExpressionResult.H"
#include "repositories/MeshesRepository.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class CommonValueExpressionDriver Declaration
\*---------------------------------------------------------------------------*/

class ExpressionDriverWriter;
class initSwakFunctionObject;
class CommonPluginFunction;

class CommonValueExpressionDriver
{
    // Private data

        //- Dictionary with all the data
        const dictionary& dict_;

        //- Empty dictionary for constructors that don't get a dictionary
        static dictionary emptyData_;

        //- Maximum depth to search for other dictionaries
        static const label maxVariableRecursionDepth_;

        //- Variables from a dictionary
        List<exprString> variableStrings_;

        //- This one is in charge for writing and restoring
        //  by using it this class doesn't have to inherit from regIOobject
        //  (which would be an overkill and inconvenient)
        autoPtr<ExpressionDriverWriter> writer_;

        //- Pointer to the "default" mesh
        static const fvMesh* defaultMeshPtr_;

        //- Getting our own pointMesh
        autoPtr<pointMesh> pMesh_;

        //- The position to which the parser progressed
        label parserLastPos_;

        //- String that allows more detailed allocation of the context
        string contextString_;

        //- String that stores the name of the other mesh
        //  (if that is going to be used)
        word otherMeshName_;

        //- List with aliases. Meant to workaround situations
        //  where fields have names that don't fit
        //  what swak accepts as a field name
        HashTable<word, word> aliases_;


protected:

    // Protected Member Functions

        //- Get the string of the context in a form that is usable for an error message
        std::string getContextString();

        const string contextString() const
        {
            return contextString_;
        }

        //- Access last position
        const label& parserLastPos() const
        {
            return parserLastPos_;
        }

        //- Access last position
        label& parserLastPos()
        {
            return parserLastPos_;
        }

        //- Get the result from another driver. Overrideable to allow mapping
        virtual autoPtr<ExpressionResult> getRemoteResult
        (
            CommonValueExpressionDriver& otherDriver
        );

        //- Is there an alias with that name
        bool hasAlias(const word& name) const;

        //- Return the actual name for an alias
        const word& getAlias(const word& name) const;

        const HashTable<word, word>& aliases() const
        {
            return aliases_;
        }

        void setAliases(const HashTable<word, word>& orig)
        {
            aliases_ = orig;
        }

        const word& resolveAlias(const word& name) const;

        const dictionary& dict() const
        {
            return dict_;
        }

        static const dictionary& emptyDict()
        {
            return emptyData_;
        }


private:

    // Private Data

        ExpressionResult result_;

        HashTable<ExpressionResult, word> variables_;

        List<StoredExpressionResult> storedVariables_;

        HashTable<DelayedExpressionResult,word> delayedVariables_;

        label specialVariablesIndex_;

        wordList globalVariableScopes_;

        HashTable<interpolationTable<scalar>, word> lines_;

        HashTable<interpolationTable<scalar>, word> lookup_;

#ifdef FOAM_HAS_INTERPOLATION2DTABLE
        HashTable<interpolation2DTable<scalar>, word> lookup2D_;
#endif

        exprString content_;

        bool trace_scanning_;

        bool trace_parsing_;

        bool variableNameIdenticalToField_;

#ifndef FOAM_DLLIBRARY_USES_STATIC_METHODS
        dlLibraryTable libraries_;
#endif


    // Private Member Functions

        void readPluginLibraries(const dictionary& dict);


protected:

    // Protected Member Functions

        const ExpressionResult& lookupGlobal(const word& name) const;

        template<class T>
        bool isGlobal
        (
            const word& name,
            bool isPoint,
            label expectedSize = -1
        ) const;

        void setGlobalScopes(const wordList &);

        //- scanner state for reentrant lexer
        void* scanner_;

        // Handling the scanner
        virtual void scan_begin() = 0;
        virtual void scan_end() = 0;

        //- update stuff that needs updating
        virtual bool update();

        //- Weights of the entries if it is not a point field
        //  parametert size is the expected size of the weight field and
        //  will be used to determine whether this is correct.
        virtual tmp<scalarField> weightsNonPoint(label size) const = 0;


public:

    // Public Member Functions

        ExpressionResult& result()
        {
            return result_;
        }

        const ExpressionResult& result() const
        {
            return result_;
        }

        const exprString& content() const
        {
            return content_;
        }

        //- Actual size of the field
        virtual label size() const = 0;

        //- Actual size of the supporting point-field
        virtual label pointSize() const = 0;

        //- Weights of the entries
        //  parametert size is the expected size of the weight field and
        //  will be used to determine whether this is correct.
        //  point determines whether this is a point field or not
        tmp<scalarField> weights(label size, bool point = false) const;

        //- Access the other mesh name
        const word& otherMeshName() const
        {
            return otherMeshName_;
        }

        word& otherMeshName()
        {
            return otherMeshName_;
        }

        //- Add an alias. Return true if the alias wasn't already there
        bool addAlias
        (
            const word& name,
            const word& field,
            bool overwrite = false
        );


protected:

    // Protected Member Functions

        //- Get the mesh we're living on
        virtual const fvMesh& mesh() const = 0;

        //- Get the mesh we're living on
        virtual const objectRegistry& db() const = 0;

        //- Get a point mesh (if needed)
        inline const pointMesh& pMesh();

        //- Get the current time
        virtual const word time() const;

        //- Set the default Mesh (if it is unset)
        static bool setUnsetDefaultMesh(const fvMesh&);

        //- Set the default Mesh
        static bool resetDefaultMesh(const fvMesh&);

        //- Get the default Mesh (if one is defined)
        static const fvMesh& getDefaultMesh();

        //- Find out which mesh is meant from the dictionary
        static const fvMesh& regionMesh
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool readIfNecessary
        );

        //- Look through the current variable values
        //  and update the stored variables
        void updateSpecialVariables(bool force = false);

        //- Read data from the writer
        virtual void getData(const dictionary&);

        //- Write data
        virtual void prepareData(dictionary&) const;

        //- Do we need a datafile written
        virtual bool hasDataToWrite() const;


private:

    // Private Data

        //- Keep fields read from disc in memory
        bool cacheReadFields_;

        //- Should we search for the field
        //  in the registry before looking on disc
        bool searchInMemory_;

        //- Should we look on disc? does not make sense in the application
        bool searchOnDisc_;

        //- Should cellSets and faceSets be cached
        //  instead of read from disc every time
        static bool cacheSets_;

        //- Use the value of the previous iteration
        //  if the old time ist asked for
        bool prevIterIsOldTime_;


protected:

    // Protected Member Functions

        inline bool cacheReadFields() const
        {
            return cacheReadFields_;
        }

        inline bool searchInMemory() const
        {
            return searchInMemory_;
        }

        inline bool searchOnDisc() const
        {
            return searchOnDisc_;
        }

        inline bool cacheSets() const
        {
            return cacheSets_;
        }

        inline bool prevIterIsOldTime() const
        {
            return prevIterIsOldTime_;
        }

        word getTypeOfFieldInternal
        (
            const fvMesh& mesh,
            const word& name
        ) const;

        word getTypeOfField(const word& name) const;

        word getTypeOfSet(const word& name) const;


public:

    // Public Member Functions

        bool isCellSet(const word& name);
        bool isCellZone(const word& name);

        bool isFaceSet(const word& name);
        bool isFaceZone(const word& name);

        bool isPointSet(const word& name);
        bool isPointZone(const word& name);


protected:

    // Protected Member Functions

        template<class T,class Mesh>
        inline tmp<T> readAndRegister
        (
            const word& name,
            const Mesh& actualMesh
        );

        enum SetOrigin
        {
            NEW,
            FILE,
            MEMORY,
            CACHE,
            INVALID
        };

        template<class T>
        inline autoPtr<T> getTopoSet
        (
            const fvMesh&,
            const word&,
            SetOrigin& origin
        ) const;

        template<class T>
        inline bool updateSet
        (
            autoPtr<T>& theSet,
            const word&,
            SetOrigin origin
        ) const;

        template<class T>
        autoPtr<T> getOrReadField
        (
            const word& name,
            bool failIfNoneFound = true,
            bool getOldTime = false
        );

        template<class T>
        autoPtr<T> getOrReadPointField
        (
            const word& name,
            bool failIfNoneFound = true,
            bool getOldTime = false
        );

        template<class T,class Mesh>
        autoPtr<T> getOrReadFieldInternal
        (
            const word& name,
            const Mesh& mesh,
            bool failIfNoneFound = true,
            bool getOldTime = false
        );

        template<class T>
        inline void correctField(T& f) const;

        typedef HashTable<int> StartupSymbols;

        template<class T>
        class SymbolTable
        :
        public StartupSymbols
        {
        public:
            SymbolTable();
        };

        virtual int startupSymbol(const word& name) = 0;

        virtual void parseInternal(int startSymbol) = 0;

        //- Only subclasses have write-access
        ExpressionResult& variable(const word& name);

        //- Expand string list to include other lists reference with #
        static List<string> expandIncludeStringList
        (
            const List<string>& orig,
            const dictionary& dict,
            const label recursionDepth
        );

        static exprString expandDictVariables
        (
            const string& orig,
            const dictionary& dict
        );


public:

    //- Runtime type information
    TypeName("CommonValueExpressionDriver");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CommonValueExpressionDriver,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict,mesh)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        CommonValueExpressionDriver,
        idName,
        (
            const word& id,
            const fvMesh& mesh,
            const objectRegistry& obr
        ),
        (id,mesh,obr)
    );

    // Constructors

        //- Construct from components
        CommonValueExpressionDriver
        (
            bool cacheReadFields = false,
            bool searchInMemory = true,
            bool searchOnDisc = false,
            const dictionary& dict = emptyData_
        );

        //- Copy construct
        CommonValueExpressionDriver(const CommonValueExpressionDriver&);

        virtual autoPtr<CommonValueExpressionDriver> clone() const = 0;

        //- Construct from a dictionary
        CommonValueExpressionDriver(const dictionary&);

        //- Return a reference to the selected value driver
        static autoPtr<CommonValueExpressionDriver> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return a reference to the selected value driver
        static autoPtr<CommonValueExpressionDriver> New
        (
            const dictionary& dict
        );

        //- Return a reference to the selected value driver
        static autoPtr<CommonValueExpressionDriver> New
        (
            const word& type,
            const word& id,
            const fvMesh& mesh,
            const objectRegistry& obr
        );


    //- Destructor
    virtual ~CommonValueExpressionDriver();


    // Member Functions

        const ExpressionResult& variable(const word& name) const;

        bool hasVariable(const word& name) const;

        void setSearchBehaviour
        (
            bool cacheReadFields,
            bool searchInMemory,
            bool searchOnDisc
        );

        template<class Type>
        tmp<Field<Type>> evaluate(const exprString& expr, bool qPoint = false);

        template<class Type>
        Type evaluateUniform(const exprString& expr, bool qPoint = false);

        void evaluateVariable(const word& name,const exprString& expr);

        void evaluateVariableRemote
        (
            const exprString& remoteExpr,
            const word& name,
            const exprString& expr
        );

        //- Allows adding uniform variables from outside
        template<class T>
        void addUniformVariable
        (
            const word& name,
            const T& val
        );

        //- Calculate the weighted average of a field
        template<class Type>
        Type calcWeightedAverage(const Field<Type>& f) const;

        //- Read information about the foreign meshes and create them
        static label readForeignMeshInfo
        (
            const dictionary& dict,
            bool verbose = true
        );

        //- Read the list of variable strings
        //  (or initialize with a single string)
        static List<exprString> readVariableStrings
        (
            const dictionary& dict,
            const word& name = "variables",
            const label recursionDepth = 0
        );

        //- Read an expression string and also do substitution
        static exprString readExpression
        (
            const word& name,
            const dictionary& dict
        );

        //- Read an expression string and also do substitution
        exprString readExpression(const word& name);

        //- Read the list of variable strings
        //  (or initialize with a single string)
        void setVariableStrings(const dictionary& dict);

        Ostream& writeVariableStrings(Ostream& out) const;

        Ostream& writeCommon(Ostream& out, bool debug = false) const;

        void addVariables(const exprString& exprList, bool clear = true);

        void addVariables
        (
            const List<exprString>& exprList,
            bool clear = true
        );

        template<class Type>
        void setResult(Field<Type>*, bool isPoint = false);

        template<class Type>
        tmp<Field<Type>> getResult(bool qPoint = false);

        void outputResult(Ostream&);

        string outputEntry();

        virtual word getResultType();

        ExpressionResult getUniform(label size, bool noWarning);

        template<template<class> class BinOp, class Type>
        inline Type getReduced
        (
            BinOp<Type> op,
            Type initial = pTraits<Type>::zero
        );

        void clearResult();

        template <class Op,class Type>
        tmp<Field<bool>> doCompare
        (
            const Field<Type>& a,
            Op op,
            const Field<Type>& b
        );

        //- This method makes sure that the interpretation
        //  of floats as logical values is consistent
        static inline bool toBool(const scalar f);

        //- Use these constants instead of literal 0 and 1
        enum BoolValues
        {
            FALSE_Value = 0,
            TRUE_Value = 1
        };

        template <class Op>
        tmp<Field<bool>> doLogicalOp
        (
            const Field<bool>& a,
            Op op,
            const Field<bool>& b
        );

        inline tmp<Field<bool>> doLogicalNot(const Field<bool>& a);

        template<class Type>
        tmp<Field<Type>> doConditional
        (
            const Field<bool>& d,
            const Field<Type>& yes,
            const Field<Type>& no
        );

        template<class Type>
        tmp<Field<Type>> makeField(Type,bool isPoint=false);

        template<class Type>
        tmp<Field<Type>> getVariable
        (
            const word& name,
            const label expectedSize
        );

        template<class Type>
        bool is(const word& name, bool isPoint = false) const;

        template<class Type>
        bool isVariable
        (
            const word& name,
            bool isPoint = false,
            label expectSize = -1
        ) const;

        template<class Type>
        bool isField(const word& name, bool isPoint = false) const;

        //- Checking for variables
        template<class T>
        bool isThere(const word& name) const;

        //- Is this a timeline
        inline bool isLine(const word& name) const;

        //- Get the value for a timeline for a specific time
        tmp<scalarField> getLine(const word& name, scalar t);

        //- Value of the timeline
        scalar getLineValue(const word& name, scalar t);

        //- Read some interpolation tables
        template<typename TableType>
        void readTables
        (
            Istream& is,
            HashTable<TableType>&,
            bool clear = true
        );

        //- Is this mesh in the repository
        bool isForeignMesh(const word& name) const;

        template<class T>
        bool isForeignField
        (
            const word& meshName,
            const word& fieldName
        ) const;

        template<class T>
        tmp<T> interpolateForeignField
        (
            const word& meshName,
            const word& fieldName,
            meshToMeshOrder theOrder =
#ifdef FOAM_NEW_MESH2MESH
            meshToMesh::imCellVolumeWeight
#else
            meshToMesh::INTERPOLATE
#endif
        ) const;

        //- Read all the interpolation tables known to man
        void readTables(const dictionary&);

        //- Read variables and tables
        //  (for objects that were not constructed by a dictionary)
        virtual void readVariablesAndTables(const dictionary&);

        //- Write some interpolation tables
        template<class TableType>
        void writeTables(Ostream& os, const HashTable<TableType>&) const;

        //- Is this a lookup table
        inline bool isLookup(const word& name) const;

        //- Get the value for a lookuptable for a set of values
        tmp<scalarField> getLookup(const word& name, const scalarField& val);

        //- Is this a 2D lookup table
        inline bool isLookup2D(const word& name) const;

        //- Get the value for a 2D lookuptable for a set of values
        tmp<scalarField> getLookup2D
        (
            const word& name,
            const scalarField& val,
            const scalarField& val2
        );

        tmp<scalarField> makeModuloField
        (
            const scalarField& a,
            const scalarField& b
        ) const;

        tmp<scalarField> makeRandomField(label seed = 0) const;

        tmp<scalarField> makeGaussRandomField(label seed = 0) const;

        vector getPositionOfMinimum
        (
            const scalarField&,
            const vectorField&
        ) const;

        vector getPositionOfMaximum
        (
            const scalarField&,
            const vectorField&
        ) const;

        const Time& runTime() const;

        tmp<vectorField> composeVectorField
        (
            const scalarField&,
            const scalarField&,
            const scalarField&
        );

        tmp<tensorField> composeTensorField
        (
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&
        );

        tmp<symmTensorField> composeSymmTensorField
        (
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&
        );

        tmp<sphericalTensorField> composeSphericalTensorField
        (
            const scalarField&
        );


        // Error handling.

        template<class Location>
        void error(const Location& l, const std::string& m);

        void error(const std::string& m);

        //- Handling the parser
        label parse(const exprString& f, const word& start = "");

        //- Handling temporary variables
        void clearVariables();

        //- Set the trace variables
        void setTrace(bool scanning, bool parsing);

        bool traceScanning() const
        {
            return trace_scanning_;
        }

        bool traceParsing() const
        {
            return trace_parsing_;
        }

        //- Create a writer for this object
        void createWriterAndRead(const word& name);

        //- Write data if apropriate
        //  Should enable exact restarts
        void tryWrite() const;

        //- Return a new plugin-function
        virtual autoPtr<CommonPluginFunction> newPluginFunction
        (
            const word& name
        ) = 0;

        //- Tests for a plugin-function
        virtual bool existsPluginFunction(const word& name) = 0;

        //- Append a string to the conext for error messages
        void appendToContext(const string& value);

        friend class ExpressionDriverWriter;

        friend class initSwakFunctionObject;

        friend class CommonPluginFunction;
        friend class exprString;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CommonValueExpressionDriver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CommonValueExpressionDriverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
