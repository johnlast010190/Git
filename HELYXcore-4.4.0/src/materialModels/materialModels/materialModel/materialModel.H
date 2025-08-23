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
    (c) 2021-2024 Engys Ltd.

Class
    Foam::materialModel

Description
    Base class for selecting material models.

SourceFiles
    materialModel.C
    materialModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef materialModel_H
#define materialModel_H

#include "materialModels/materialTables/materialTables.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class materialModel Declaration
\*---------------------------------------------------------------------------*/

class materialModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Materials object registry
        const objectRegistry& matObr_;

        //- Material models tables
        materialTables& materialTables_;

        //- Phase name
        const word phaseName_;

        //- Specie name
        const word specieName_;

        //- Model name (in object registry)
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Object registry reference
        const objectRegistry& obr_;

        //- Pointer to dictionary/sub-dictionary
        const dictionary* dict_;

        //- Dependance ptr list
        depList dep_;

        //- Pointer list to the models
        matScalarUTable sMod_;

        //- Pointer list to the models
        matVectorUTable vMod_;

        //- Pointer list to the models
        matTensorUTable tMod_;

        //- Pointer to the temperature offset by reference
        const refScalarField* T_;

        //- Pointer to the pressure offset by the reference
        const refScalarField* p_;


    // Protected member functions

        //- Look up an optional dependent model and return if it is constant
        template<class Type>
        bool isModelConst() const;

        //- Initialise or construct and create pointer to requested field
        void initialise(const word& fieldName)
        {
            if (fieldName == "p")
            {
                p_ = constructOrReturnRefFieldPtr<scalar>("p");
            }
            else if (fieldName == "T")
            {
                T_ = constructOrReturnRefFieldPtr<scalar>("T");
            }
            else
            {
                FatalErrorInFunction
                    << "Requested field " << fieldName << " is not supported"
                    << exit(FatalError);
            }
        }


public:

    //- Runtime type information
    TypeName("materialModel");

    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        materialModel,
        dictionary,
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        ),
        (obr, dict, phaseName, specieName, name)
    );


    // Selectors

        //- Return reference to selected rho model
        static autoPtr<materialModel> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& modelType,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Constructors

        //- Construct from components
        materialModel
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name = word::null
        );


    //- Destructor
    virtual ~materialModel();


    // Public member functions

        //- Constructs phase filed name
        word phasePropertyName
        (
            const word fieldName,
            const word& phaseName = word::null
        ) const
        {
            return materialTables_.phasePropertyName(fieldName, phaseName);
        }

        //- Helper function to fill depenedency list with multiple models
        //  All of them are set as compulsory.
        void fill
        (
            const word& checkModel,
            const word& modelName,
            const wordList& requiredModels,
            const labelList& enumInds,
            const label depNumber
        );

        //- Helper function to fill dependency list with a single model
        void fill
        (
            const word& checkModel,
            const word& modelName,
            const word& requiredModel,
            const label enumInd,
            const label depNumber
        );

        //- Fill dependencies
        //  First check is if checkModel == modelName
        //  Than dependencies from modelName are added
        void fill
        (
            const word& checkModel,         // Check if the model name corresponds to this
            const word& modelName,          // Model for which dependencies are filled
            const wordList& requiredModels, // Required models for modelName
            const labelList& enumInds,     // Indexes to put dependent models in
            const boolList& compulsory,     // Each model has a flag if it is required or not
            const label depNumber           // Slot number for overall dep_ list
        );

        //- Helper function to fill dependency list with a single model
        void fill
        (
            const word& checkModel,
            const wordList& modelNames,
            const wordList& requiredModels,
            const labelList& enumInds
        );

        //- Fill dependencies
        void fill
        (
            const word& checkModel,
            const wordList& modList,
            const List<wordList>& dependencies,
            const List<labelList>& enumInds
        );

        void fill
        (
            const word& checkModel,
            const wordList& modList,
            const wordList& dependencies,
            const labelList& enumInds,
            const boolList& compulsory
        );

        void fill
        (
            const word& checkModel,
            const wordList& modList,
            const List<wordList>& dependencies,
            const List<labelList>& enumInds,
            const List<boolList>& compulsory
        );

        //- Constructs phase filed name
        void updateDictPtr();

        const dictionary& materialsDict() const;

        //- Get field pointer with the nullptr check
        template<class Type>
        VolField<Type>* lookupPtr
        (
            const word& fieldSubName,
            const word& phaseName = word::null
        );

        //- Constructiong or returning pointer to reference field
        template<class Type>
        const referenceFields<Type>* constructOrReturnRefFieldPtr
        (
            const word& fieldName
        ) const;

        //- Mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Update table
        virtual void updateTable(const word&) = 0;

        //- Type independent pointer list
        virtual const depList& dep() const
        {
            return dep_;
        }

        //- Casting of scalar models
        virtual baseModels<scalar>* castScalarModel(const word&)
        {
            return nullptr;
        }

        //- Casting of vector models
        virtual baseModels<vector>* castVectorModel(const word&)
        {
            return nullptr;
        }

        //- Casting of tensor models
        virtual baseModels<tensor>* castTensorModel(const word&)
        {
            return nullptr;
        }

        //- Explicitly update scalar field
        virtual bool updateScalarField
        (
            const word& fieldName,
            const volScalarField& volField
        )
        {
            return false;
        };

        //- Explicitly update vector field
        virtual bool updateVectorField
        (
            const word& fieldName,
            const volVectorField& volField
        )
        {
            return false;
        };

        //- Explicitly update tensor field
        virtual bool updateTensorField
        (
            const word& fieldName,
            const volTensorField& volField
        )
        {
            return false;
        };

        //- Is it mixture model
        virtual bool mixture() const
        {
            return false;
        }

        //- Write for data
        virtual bool writeData(Ostream& os) const override
        {
            NotImplemented;
        }

        //- Write
        virtual bool write(const bool write) const override
        {
            NotImplemented;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "materialModelTemplates.C"
#endif

#endif


// ************************************************************************* //
