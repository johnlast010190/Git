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
    (c) 2016-2024 Engys Ltd.
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::radiationModel

Description
    Top level model for radiation modelling

SourceFiles
    radiationModel.C
    radiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModel_H
#define radiationModel_H

#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatricesFwd.H"
#include "primitives/bools/Switch/Switch.H"
#include "transmissivityModels/transmissivityModel/transmissivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;
class fvMesh;

namespace radiationModels
{

// Forward declaration of classes
class absorptionEmissionModel;
class scatterModel;
class sootModel;
}

/*---------------------------------------------------------------------------*\
                       Class radiationModel Declaration
\*---------------------------------------------------------------------------*/

class radiationModel
:
    public IOdictionary
{
public:

    // Static data

        //- Name of dictionary
        static const word dictName;

        //- Static name external radiative fluxes
        static const word externalRadHeatFieldName_;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Offset for absolute temperature
        const dimensionedScalar TRef_;

        //- Participating model on/off flag
        Switch participating_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<radiationModels::absorptionEmissionModel>
                absorptionEmission_;

            //- Scatter model
            autoPtr<radiationModels::scatterModel> scatter_;

            //- Soot model
            autoPtr<radiationModels::sootModel> soot_;

            //- Transmissivity model
            autoPtr<radiationModels::transmissivityModel> transmissivity_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const objectRegistry& obr) const;

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("radiationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            radiationModel,
            T,
            (
                const volScalarField& T,
                const dimensionedScalar& TRef
            ),
            (T, TRef)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            radiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T,
                const dimensionedScalar& TRef
            ),
            (dict, T, TRef)
        );


    // Constructors

        //- Null constructor
        radiationModel(const volScalarField& T, const dimensionedScalar& TRef);

        //- Construct from components
        radiationModel
        (
            const word& type,
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Construct from components
        radiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Disallow default bitwise copy construction
        radiationModel(const radiationModel&) = delete;


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<radiationModel> New
        (
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Return a reference to the selected radiation model
        static autoPtr<radiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Lookup or create on the object registry
        static radiationModel& lookupOrCreate
        (
            const volScalarField& T,
            const dimensionedScalar& TRef
        );


    //- Destructor
    virtual ~radiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine. Return whether recalculation
            //  was done (i.e. whether calculate was called)
            virtual bool correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Participating model on/off flag
            const Switch participating() const
            {
                return participating_;
            }

            //- Return access to the temperature field
            const volScalarField& T() const noexcept
            {
                return T_;
            }

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const basicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

//            virtual label nBands() const = 0;

            //- Access to absorptionEmission model
            const radiationModels::absorptionEmissionModel&
            absorptionEmission() const;

            //- Access to transmissivity Model
            const radiationModels::transmissivityModel& transmissivity() const;

            //- Access to soot Model
            const radiationModels::sootModel& soot() const;

            //- Return temperature offset
            const dimensionedScalar& TRef() const
            {
                return TRef_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        radiationModel,                                                        \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        radiationModel,                                                        \
        model,                                                                 \
        T                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
