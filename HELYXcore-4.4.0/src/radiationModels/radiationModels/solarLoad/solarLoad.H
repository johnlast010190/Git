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
    (c) 2015-2019 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::radiationModels::solarLoad

Description

    The solar load radiation model includes Sun primary hits, their
    reflective fluxes and diffusive sky radiative fluxes.

    The primary hit rays are calculated using a face shading algorithm.
    The reflected fluxes are considered diffusive and use a view factors method
    to deposit the energy on "visible" walls. The sky diffusive radiation for
    horizontal and vertical walls is calculated following the Fair Weather
    Conditions Method from the ASHRAE Handbook.

    By default the energy is included in cells adjacent to the patches into
    the energy Equation (wallCoupled = false). On coupled patches the flux is
    by default added to the wall and considered into the solid
    (solidCoupled = true).

    The reflected fluxes uses a grey absoprtion/emission model wich is weighted
    by the spectral distribution. The flag useVFbeamToDiffuse should be
    switched on and the view factors should be calculated using the
    'viewFactorsGen' application.

    The solarLoad model can be used in conjuntion with fvDOM and viewFactor
    radiation models but only using a single band spectrum. On the
    corresponding BC's for these models the flag "solarLoad" must be set to
    true.


SourceFiles
    solarLoad.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelsolarLoad_H
#define radiationModelsolarLoad_H

#include "radiationModels/radiationModel/radiationModel.H"
#include "fvMesh/singleCellFvMesh/singleCellFvMesh.H"
#include "primitives/Scalar/lists/scalarListIOList.H"
#include "fields/volFields/volFields.H"
#include "radiationModels/solarLoad/faceShading/faceShading.H"
#include "solarCalculator/solarCalculator.H"
#include "meshes/polyMesh/polyDistributionMap/IOdistributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class solarLoad Declaration
\*---------------------------------------------------------------------------*/

class solarLoad
:
    public radiationModel
{
protected:

    // Static data

        //- Static name for view factor walls
        static const word viewFactorWalls;


private:

    // Private data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Coarse mesh
        autoPtr<singleCellFvMesh> coarseMesh_;

        //- Net radiative heat flux [W/m2]
        volScalarField qr_;

        //- Secondary solar radiative heat flux [W/m2]
        volScalarField qsecondRad_;

        //- Direct hit faces Ids
        autoPtr<faceShading> hitFaces_;

        //- Constant source term
        DimensionedField<scalar, volMesh> Ru_;

        //- Solar calculator
        solarCalculator solarCalc_;

        //- Vertical direction (Default is g vector)
        vector verticalDir_;

        //- Include diffuse reflected heat fluxes from direct heat flux
        bool useVFbeamToDiffuse_;

        //- Selected patches to apply solar load
        labelList includePatches_;

        //- Chached coarse to fine mapping for coarse mesh
        List<labelListList> coarseToFine_;

        //-Number of bands
        label nBands_;

        //- Spectral distribution for the integrated solar heat flux
        scalarList spectralDistribution_;

        //- Time-dependent spectral distributions
        autoPtr<Function1<scalarField>> spectralDistributions_;

        //- Primary solar radiative heat flux per band [W/m2]
        PtrList<volScalarField> qprimaryRad_;

        //- Map distribute
        autoPtr<IOdistributionMap> map_;

        //- Face-compact map
        labelListIOList visibleFaceFaces_;

        //- Couple solids through mapped boundary patch using qr (default:true)
        bool solidCoupled_;

        //- Couple wall patches using qr (default:false)
        bool wallCoupled_;

        //- Absorptivity list
        List<List<tmp<scalarField>>> absorptivity_;

        //- Update absorptivity
        bool updateAbsorptivity_;

        //- First iteration
        bool firstIter_;

        //- Update Sun position index
        label updateTimeIndex_;


    // Private Member Functions

        //- Initialise
        void initialise(const dictionary&);

        //- Update direct hit faces radiation
        void updateDirectHitRadiation(const labelList&, const labelHashSet&);

        //- Calculate diffusive heat flux
        void calculateQdiff(const labelHashSet&, const labelHashSet&);

        //- Update Sky diffusive radiation
        void updateSkyDiffusiveRadiation
        (
            const labelHashSet&,
            const labelHashSet&
        );

        //- Update hit faces
        bool updateHitFaces();

        //- Update absorptivity
        void updateAbsorptivity(const labelHashSet& includePatches);


public:

    //- Runtime type information
    TypeName("solarLoad");


    // Constructors

        //- Construct from volScalarField
        solarLoad(const volScalarField& T, const dimensionedScalar& TRef);

        //- Construct from dictionary and volScalarField
        solarLoad
        (
            const dictionary& dict,
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Constructor from local components. Does not create a radiationModel.
        //  radWallFieldName is the solar heat field name
        solarLoad
        (
            const dictionary& dict,
            const volScalarField& T,
            const dimensionedScalar& TRef,
            const word radWallFieldName
        );

        //- Disallow default bitwise copy construction
        solarLoad(const solarLoad&) = delete;


    //- Destructor
    virtual ~solarLoad();


    // Member functions

        //- Solve
        void calculate();

        //- Read radiation properties dictionary
        bool read();

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        virtual tmp<DimensionedField<scalar, volMesh>> Ru() const;

        //- Number of bands
        label nBands() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solarLoad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
