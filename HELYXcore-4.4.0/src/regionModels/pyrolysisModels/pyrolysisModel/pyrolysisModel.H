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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::regionModels::pyrolysisModel

Description
    Base class for pyrolysis models

SourceFiles
    pyrolysisModelI.H
    pyrolysisModel.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisModel_H
#define pyrolysisModel_H

#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "fields/volFields/volFieldsFwd.H"
#include "regionModel/regionModel/regionModel.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/ints/lists/labelList.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "multiphaseThermo/multiphaseThermo.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class Time;

namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class pyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisModel
:
    public regionModel
{
    // Private Member Functions

        //- Initialise the region
        void initialise(const word& modelType);


protected:

    // Protected data

        //- Pimple solution control
        pimpleControl pimple_;

        //- Pointer to multiphase thermo
        multiphaseThermo *thermoPtr_;

        //- Pointer to the solid chemistry model
        mutable solidThermo* solidThermo_;

        //- Reference to the solid chemistry model
        mutable autoPtr<basicSolidChemistryModel> solidChemistry_;

        //- Pointer to radiation model
        autoPtr<radiationModel> radiation_;

        //- Global face IDs
        labelListList boundaryFaceFaces_;

        //- Global cell IDs
        labelListList boundaryFaceCells_;

        //- Global boundary face IDs oppossite coupled patch
        labelList boundaryFaceOppositeFace_;

        //- Number of layers in the region
        label nLayers_;

        //- Face area magnitude normal to patch
        autoPtr<surfaceScalarField> nMagSfPtr_;

        //- Flag to allow mesh movement
        Switch moveMesh_;

        //- Maximum diffussivity
        scalar maxDiff_;


    // Protected Member Functions

        //- Read control parameters
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);

        //- Return the global face IDs
        const labelListList& boundaryFaceFaces() const
        {
            return boundaryFaceFaces_;
        }

        //- Return the global cell IDs
        const labelListList& boundaryFaceCells() const
        {
            return boundaryFaceCells_;
        }

        //- Return the global boundary face IDs oppossite coupled patch
        const labelList& boundaryFaceOppositeFace() const
        {
            return boundaryFaceOppositeFace_;
        }

        //- Return the face area magnitudes / [m2]
        const surfaceScalarField& nMagSf() const
        {
            return nMagSfPtr_();
        }

        //- Return the number of layers in the region
        inline label nLayers() const
        {
            return nLayers_;
        }

        //- Move mesh points according to change in cell volumes
        //  Returns map ordered by cell where 1 = cell moved, 0 = cell unchanged
        tmp<labelField> moveMesh
        (
            const scalarList& deltaV,
            const scalar minDelta = 0.0
        );


public:

    //- Runtime type information
    TypeName("pyrolysisModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pyrolysisModel,
            mesh,
            (
                const word& modelType,
                const fvMesh& mesh,
                const word& regionType
            ),
            (modelType, mesh, regionType)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pyrolysisModel,
            dictionary,
            (
                const word& modelType,
                const fvMesh& mesh,
                const dictionary& dict,
                const word& regionType
            ),
            (modelType, mesh, dict, regionType)
        );


    // Constructors

        //- Construct from type name and mesh
        pyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name and mesh and dictionary
        pyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );

        //- Disallow default bitwise copy construction
        pyrolysisModel(const pyrolysisModel&) = delete;

        //- Return clone
        autoPtr<pyrolysisModel> clone() const
        {
            NotImplemented;
        }


    // Selectors

        //- Return a reference to the selected pyrolysis model
        static autoPtr<pyrolysisModel> New
        (
            const fvMesh& mesh,
            const word& regionType = "pyrolysis"
        );

        //- Return a reference to a named selected pyrolysis model
        static autoPtr<pyrolysisModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType = "pyrolysis"
        );


    //- Destructor
    virtual ~pyrolysisModel();


    // Member Functions

        //- Return the total gas mass flux to primary region [kg/m2/s]
        virtual const surfaceScalarField& phiGas() const
        {
            NotImplemented;
        }

        //- External hook to add mass to the primary region
        virtual scalar addMassSources(const label patchi, const label facei);

        //- Mean diffusion number of the solid region
        virtual scalar solidRegionDiffNo() const;

        //- Return max diffusivity allowed in the solid
        virtual scalar maxDiff() const;

        //- Return reference to solid thermo
        virtual const solidThermo& sThermo() const
        {
            return *solidThermo_;
        }

        //- Density
        virtual volScalarField& rho()
        {
            return (*solidThermo_).rho();
        }

        //- Density
        virtual const tmp<volScalarField> rho() const
        {
            return (*solidThermo_).rho();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pyrolysisModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
