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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2024 Engys Ltd

Class
    Foam::radiationModels::radiativeIntensityRay

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    radiativeIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeIntensityRay_H
#define radiativeIntensityRay_H

#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "radiationModels/fvDOM/blackBodyEmission/blackBodyEmission.H"
#include "finiteVolume/fvm/fvm.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

// Forward declaration of classes
class fvDOM;

/*---------------------------------------------------------------------------*\
                    Class radiativeIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class radiativeIntensityRay
{
public:

    static const word intensityPrefix;


private:

    // Private data

        //- Refence to the owner fvDOM object
        const fvDOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Absorption/emission model
        const absorptionEmissionModel& absorptionEmission_;

        //- Black body
        const blackBodyEmission& blackBody_;

        //- Total radiative intensity / [W/m2]
         volScalarField I_;

        //- Total (net) radiative heat flux on boundary
        autoPtr<volScalarField::Boundary> qr_;

        //- Incident (absorbed) radiative heat flux on boundary
        autoPtr<volScalarField::Boundary> qin_;

        //- Emitted radiative heat flux on boundary
        autoPtr<volScalarField::Boundary> qem_;

        //- Irradiation heat flux on boundary
        autoPtr<volScalarField::Boundary> qg_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to radiative intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;

        //- solution order for non-participating radiation solver
        autoPtr<labelList> solutionOrderPtr_;

        //- Global ray id - incremented in constructor
        static label rayId;

        //- My ray Id
        label myRayId_;


public:

    // Constructors

        //- Construct form components
        radiativeIntensityRay
        (
            const fvDOM& dom,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const absorptionEmissionModel& absEmmModel_,
            const blackBodyEmission& blackBody,
            const label rayId
        );

        //- Disallow default bitwise copy construction
        radiativeIntensityRay(const radiativeIntensityRay&) = delete;


    //- Destructor
    ~radiativeIntensityRay();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

            //- Initialise the ray in i direction
            void init
            (
                const scalar phi,
                const scalar theta,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add radiative intensities from all the bands
            void addIntensity();

            //- calculate patch heat flux due to emitted radiation
            tmp<scalarField> emittedRadiantIntensity
            (
                label patchI,
                const scalarField& TpAbs
            ) const;


        // Access

            //- Return intensity
            inline const volScalarField& I() const;

            //- Return cell integrated Qr
            tmp<volScalarField> cellQr() const;

            //- Return const access to the boundary total (net) heat flux
            inline const volScalarField::Boundary& qr() const;

            //- Return non-const access to the boundary total (net) heat flux
            inline volScalarField::Boundary& qr();

            //- Return non-const access to the boundary incident (absorbed) heat flux
            inline volScalarField::Boundary& qin();

            //- Return const access to the boundary incident (absorbed) heat flux
            inline const volScalarField::Boundary& qin() const;

            //- Return non-const access to the boundary emitted heat flux
            inline volScalarField::Boundary& qem();

            //- Return const access to the boundary emmited heat flux
            inline const volScalarField::Boundary& qem() const;

            //- Return non-const access to the boundary irradiation heat flux
            inline volScalarField::Boundary& qg();

            //- Return const access to the boundary irradiation heat flux
            inline const volScalarField::Boundary& qg() const;

            //- Return direction
            inline const vector& d() const;

            inline vector& d();

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            inline vector& dAve();

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the name of the ray
            inline word rayName() const;

            //- Return the radiative intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiativeIntensityRay&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiationModels/fvDOM/radiativeIntensityRay/radiativeIntensityRayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
