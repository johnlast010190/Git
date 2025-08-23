/*---------------------------------------------------------------------------* \
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
    (c) 2023 Engys Ltd.

Class
    Foam::equationRelax

Description
    helper class to relax equations

Usage

    Example of the equation relax class usage from the solverObject
    \verbatim
    <solverObjectName>
    {
        ...
        localRelax true;
        localRelaxCoeffs // if localRelax is true
        {
            //- if dict is empty, default values are used
            relaxType Courant;
            minCourant 0.1;
            maxCourant 1e5;
            minCourantRelax 0.99999;
            maxCourantRelax 0.9;
        }

        meshQualityRelax 1;
        meshQualityControls // if meshQualityRelax is less than 1
        {
            //- if dict is empty, default values are used
            nonOrthoThreshold 60;
            skewnessThreshold 0.9;
        }
        procRelax 1; // relaxes all processor patchCells
        localCellProcRelax 1; // relaxes cells that are on walls and have less
                              // than two internal faces. This is usually
                              // happens when the mesh is decomposed in a not
                              // optimal way
        ...
    }
    \endverbatim

SourceFiles
    equationRelax.C

\*---------------------------------------------------------------------------*/
#ifndef equationRelax_H
#define equationRelax_H

#include "cfdTools/general/include/fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class equationRelax Declaration
\*---------------------------------------------------------------------------*/

class equationRelax
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the dictionary
        const dictionary& dict_;

        //- Time mode -- pseudoTime (steady) or physical (transient)
        const word timeMode_;

        //- local relaxation Switch
        const Switch localRelax_;

        //- PseudoTime - relax type
        word relaxType_;

        //- Coeffs needed for Courant based relax type
        //- The user specifies a range of Courant [minC, maxC]
        //- and a relax that he wants to specify there [minCRelax, maxCRelax]
        //  The maxCourantRelax should value should be lower than the
        //  minCourantRelax value since we want more relax to be applied in
        //  higher Courant values
        autoPtr<Function1<scalar>> minC_;
        autoPtr<Function1<scalar>> maxC_;
        autoPtr<Function1<scalar>> minCRelax_;
        autoPtr<Function1<scalar>> maxCRelax_;

        //- Relax processor patch cells
        const scalar procRelax_;

        //- Relax wall patch cells
        const scalar wallRelax_;

        //- Relax based on mesh quality
            //- Relax bad quality cells
            const scalar meshQualityRelax_;

            //- nonOrthogonality threshold
            scalar nonOrthoThreshold_;

            //- skewness threshold
            scalar skewnessThreshold_;

            //- bad quality cell marker
            labelList badCellQualityMarker_;

        //- Relax on bad decomposition cells
        const scalar localCellProcRelax_;

        //- bad decomposition cell marker
        labelList badDecompositionCellMarker_;

        //- Relax cells with AMI and processor faces
        const scalar localCellProcAMIRelax_;

        //- Cells that have AMI and processor faces
        labelList procAndAMICellMarker_;

        //- Calculate the point cell delta;
        autoPtr<vectorField> pcdelta_;


    // Private Functions

        //- Initialize
        void initialize();

        //- Set Courant relax coeffs
        void setCourantRelaxCoeffs();

        //- Use non-uniform momentum relaxation
        void nonUniformMomentumRelax(fvBlockMatrix<vector>&) const;

        //- Compute local relaxation factor for momentum based on relaxation type
        tmp<vectorField> computeLocalRelaxFactor
        (
            const volVectorField&
        ) const;

        //- Compute Courant No for each dimension
        //  Co_i = U_i*dt/dx_i
        tmp<vectorField> computeCourantNo(const volVectorField&) const;

        //-assign deltas
        void calculateCellDeltas();

        //- Mark bad quality cells
        labelList badCellQualityMarker(bool printInfo=true) const;

        //- Mark cells with placed in multiple boundaries
        labelList badDecompositionCellMarker(bool printInfo=true) const;

        //- Mark cells with AMI and processor faces
        labelList cellsWithProcAndAMIFaces(bool printInfo=true) const;

public:

// Constructors

    // - Construct from components
    equationRelax
    (
        const fvMesh& mesh,
        const dictionary&,
        const bool=true
    );

    //- Destructor
    virtual ~equationRelax() = default;


    // Public member functions

        //- Access
        bool isTransient() const
        {
            return timeMode_ == "physical";
        }

        //- Momentum relax
        void momentumRelax(fvBlockMatrix<vector>&) const;

        //- Update bad cells when mesh is changed
        void updateCellMarkers();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
