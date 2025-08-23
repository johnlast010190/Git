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
    (c) 2015 OpenCFD Ltd.
    (c) 2010-2016 Engys Ltd.

Class
    Foam::boundaryRadiationPropertiesPatch

Description
    Common functions to emissivity, transmissivity and transmission source.
    Emissivity and transmissivity supplied from lookup into a
    dictionary or calculated by the solidThermo:

    - 'lookup' : lookup volScalarField with name
    - 'solidRadiation' : use solidThermo
    - 'model' : use a model

    Transmission source currently hardcoded as time dependent ambient
    temperature or heat flux

    - Ta [K]
    - Qa [W/m2]

SourceFiles
    boundaryRadiationPropertiesPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationPropertiesPatch_H
#define boundaryRadiationPropertiesPatch_H

#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/enums/NamedEnum.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "fields/fvPatchFields/basic/calculated/calculatedFvPatchFields.H"
#include "transmissivityModels/transmissivityModel/transmissivityModel.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "primitives/functions/Function1/Function1/Function1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
/*---------------------------------------------------------------------------*\
             Class boundaryRadiationPropertiesPatch Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationPropertiesPatch
{
public:

    // - Type of method
    enum methodType
    {
        SOLIDRADIATION,
        LOOKUP,
        MODEL
    };

    enum modelType
    {
        TRANSPARENT,
        OPAQUE
    };


private:

    // Private data

        static const NamedEnum<methodType, 3> methodTypeNames_;

        //- How to get property
        const methodType method_;

        static const NamedEnum<modelType, 2> modelTypeNames_;

        //- How to get property
        const modelType model_;

        //- Dictionary
        const dictionary dict_;

        //- Absorption/emission model
        autoPtr<absorptionEmissionModel> absorptionEmission_;

        //- transmissivity model
        autoPtr<transmissivityModel> transmissivity_;

        //- transmission source model
        // autoPtr<transmissionSource> transmissionSource_;

        //- Time dependent Function1 for ambient temperature
        mutable autoPtr<Function1<scalar>> Ta_;

        //- Time dependent Function1 for external rad. heat
        mutable autoPtr<Function1<scalar>> Qa_;

        //- reference to fvPatch
        const polyPatch& patch_;

        //- Temperature offset
        const scalar TRef_;


    // Private functions

        //- Return nbr patch index
        label nbrPatchIndex() const;

        //- Return nbr mesh
        const fvMesh& nbrRegion() const;


public:


    // Constructors


        //- Construct from patch, internal field and dictionary
        boundaryRadiationPropertiesPatch
        (
            const polyPatch&,
            const dictionary&,
            const scalar& TRef
        );

         //- Destructor
        virtual ~boundaryRadiationPropertiesPatch();


    // Member functions

        //- Method to obtain properties
        word method() const
        {
            return methodTypeNames_[method_];
        }


        //- Return absorptionEmissionModel
        const absorptionEmissionModel& absorptionEmission() const;

        //- Return transmissivityModel
        const transmissivityModel& transmissiveModel() const;

        //- Calculate corresponding emissivity field for bandI
        tmp<scalarField> emissivity(const label bandI = 0) const;

        //- Calculate corresponding absorptivity field for bandI
        tmp<scalarField> absorptivity(const label bandI = 0) const;

        //- Calculate corresponding transmissivity field for bandI
        tmp<scalarField> transmissivity(const label bandI = 0, bool useSolarTransmissivity = false) const;

        //- Calculate corresponding reflectivity field
        tmp<scalarField> reflectivity(const label bandI = 0) const;

        //- Calculate emitted radiant power for given temperature
        tmp<scalarField> emittedRadiantFlux
        (
            const scalarField& T, const label bandI = 0
        ) const;

        //- Calculate transmitted radiant power from environmental temperature
        //*** This should be generalised in future to allow different kinds
        // and combinations of transmission sources
        tmp<scalarField> radiantTransmissionSource
        (
            const label bandI = 0
        ) const;

        //- Return scalar ambient temperature
        scalar Tambient() const;

        //- Return scalar ambient rad. heat
        scalar Qambient() const;


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiationModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
