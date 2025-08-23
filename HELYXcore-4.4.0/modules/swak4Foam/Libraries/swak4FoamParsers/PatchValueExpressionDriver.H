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
    Foam::PatchValueExpressionDriver

Description

SourceFiles
    PatchValueExpressionDriverI.H
    PatchValueExpressionDriver.C

Contributors/Copyright:
    2009-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef PatchValueExpressionDriver_H
#define PatchValueExpressionDriver_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ExpressionResult/ExpressionResult.H"

#include "CommonValueExpressionDriver.H"

#include "interpolation/interpolation/interpolation/interpolation.H"

namespace Foam
{
    class PatchValueExpressionDriver;
}

using Foam::PatchValueExpressionDriver;

#include "PatchValueExpressionParser.tab.hh"

namespace parserPatch
{
    class location;
    class PatchValueExpressionParser;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchValueExpressionDriver Declaration
\*---------------------------------------------------------------------------*/

class PatchValueExpressionDriver
:
    public CommonValueExpressionDriver
{
    // Private data

        //- Registry
        const objectRegistry& obr_;

        //- The patch we're living on
        const fvPatch& patch_;

        //- Dictionary with the interpolation-schemes for mapping
        dictionary mappingInterpolationSchemes_;

        static word driverName_;

        //- HashPtrTable of interpolations
        HashPtrTable<interpolation<scalar>> interpolationScalar_;
        HashPtrTable<interpolation<vector>> interpolationVector_;
        HashPtrTable<interpolation<tensor>> interpolationTensor_;
        HashPtrTable<interpolation<symmTensor>> interpolationSymmTensor_;
        HashPtrTable<interpolation<sphericalTensor>>
            interpolationSphericalTensor_;


    // Private Member Functions

        //- Give back an interpolation-scheme specified by name and type
        template<class Type>
        const interpolation<Type>& getInterpolation(const VolField<Type>& psi);

        //- Get the specification of the interpolations by name
        const word getInterpolationScheme(const word& name);

        //- Table with the interpolation scheme (depending on the type)
        template <class Type>
        HashPtrTable<interpolation<Type>>& interpolations();

        //- Map a field that is as long as the patch
        template<class T>
        tmp<Field<T>> mapField(tmp<Field<T>> orig);


protected:

    // Protected Member Functions

        //- Handling the scanner
        void scan_begin();
        void scan_end();

        //- Set the state of the lexer to scanning vector components
        void startVectorComponent();

        //- Set the state of the lexer to scanning tensor components
        void startTensorComponent();

        //- Set the state to "eat characters"
        void startEatCharacters();

        //- Handling the parser
        void parseInternal(int startToken);

        //- Execute a plugin-function
        template<class T>
        autoPtr<Field<T>> evaluatePluginFunction
        (
            const word& name,
            const parserPatch::location& loc,
            int& scanned,
            bool isPoint = false
        );

        static const SymbolTable<PatchValueExpressionDriver>& symbolTable();

        int startupSymbol(const word& name);

        //- Check whether the patch maps from another one
        autoPtr<ExpressionResult> getRemoteResult
        (
            CommonValueExpressionDriver& otherDriver
        );

        //- Weights of the entries if it is not a point field
        //  Parametert size is the expected size of the weight field and
        //  will be used to determine whether this is correct.
        virtual tmp<scalarField> weightsNonPoint(label size) const;


public:

    //- Runtime type information
    TypeName("PatchValueExpressionDriver");

    //- Entity this driver operates on
    typedef fvPatch EntityType;

    static const word& driverName()
    {
        return driverName_;
    }


    // Constructors

        //- Construct from components
        PatchValueExpressionDriver(const fvPatch& patch_);

        //- Construct from components
        PatchValueExpressionDriver
        (
            const dictionary& dict,
            const fvPatch& patch_,
            const objectRegistry& obr
        );

        //- Construct from a dictionary
        PatchValueExpressionDriver(const dictionary& dict, const fvMesh& mesh);

        //- Construct from and an old Version
        PatchValueExpressionDriver
        (
            const fvPatch& patch_,
            const PatchValueExpressionDriver& old
        );

        //- Construct from an id
        PatchValueExpressionDriver
        (
            const word& id,
            const fvMesh& mesh,
            const objectRegistry& obr
        );

        //- Copy construct
        PatchValueExpressionDriver(const PatchValueExpressionDriver&);

        //- Clone
        autoPtr<CommonValueExpressionDriver> clone() const
        {
            return autoPtr<CommonValueExpressionDriver>
            (
                new PatchValueExpressionDriver(*this)
            );
        }

    //- Destructor
    ~PatchValueExpressionDriver();


    // Member Functions

        //- Return the patch
        const fvPatch& patch() const
        {
            return patch_;
        }

        //- Actual size of the field
        virtual label size() const;

        //- Actual size of the supporting point-field
        virtual label pointSize() const;

        //- Get the mesh we're living on
        virtual const fvMesh& mesh() const;

        //- Get registry
        virtual const objectRegistry& db() const;

        tmp<vectorField> makePositionField();

        tmp<vectorField> makePointField();

        tmp<vectorField> makeFaceNormalField();

        tmp<vectorField> makeFaceAreaField();

        tmp<vectorField> makeCellNeighbourField();

        tmp<vectorField> makeDeltaField();

        tmp<scalarField> makeWeightsField();

        tmp<scalarField> makeFaceIdField();

        tmp<scalarField> makeNearDistField();

        template<class Type>
        tmp<Field<Type>> getField(const word& name);

        template<class Type>
        tmp<Field<Type>> getOldTimeField(const word& name);

        template<class Type>
        tmp<Field<Type>> getSurfaceNormalField(const word& name);

        template<class Type>
        tmp<Field<Type>> getPatchInternalField(const word& name);

        template<class Type>
        tmp<Field<Type>> getPatchNeighbourField(const word& name);

        template<class Type>
        tmp<Field<Type>> toFace(const Field<Type>&);

        template<class Type>
        tmp<Field<Type>> toPoint(const Field<Type>&);

        //- Return a new plugin-function
        virtual autoPtr<CommonPluginFunction> newPluginFunction
        (
            const word& name
        );

        //- Tests for a plugin-function
        virtual bool existsPluginFunction(const word& name);

        //- For a mapped patch returns the mapped field from "the other side""
        template<class T>
        tmp<Field<T>> mapField(const word& name, bool internalField = false);

        //- For a mapped point-patch returns
        //  the mapped field from "the other side""
        template<class T>
        tmp<Field<T>> mapPointField(const word& name);


    //- Allow access of protected methods
    friend class parserPatch::PatchValueExpressionParser;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PatchValueExpressionDriver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchValueExpressionDriverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
