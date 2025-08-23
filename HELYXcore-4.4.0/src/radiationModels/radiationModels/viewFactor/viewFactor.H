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
    (c) 2024 Engys Ltd.

Class
    Foam::radiationModels::viewFactor

Description
    View factor radiation model. The system solved is: C q = b
    where:
            Cij  = deltaij/Ej - (1/Ej - 1)Fij
            q    = heat flux
            b    = A eb - Ho
    and:
            eb   = sigma*T^4
            Ej   = emissivity
            Aij  = deltaij - Fij
            Fij  = view factor matrix


SourceFiles
    viewFactor.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelviewFactor_H
#define radiationModelviewFactor_H

#include "radiationModels/radiationModel/radiationModel.H"
#include "fvMesh/singleCellFvMesh/singleCellFvMesh.H"
#include "matrices/scalarMatrices/scalarMatrices.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "primitives/Scalar/lists/scalarListIOList.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/polyDistributionMap/IOdistributionMap.H"
#include "radiationModels/solarLoad/solarLoad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                         Class viewFactor Declaration
\*---------------------------------------------------------------------------*/

class viewFactor
:
    public radiationModel
{
protected:

// Static data

    //- Static name for view factor walls
    static const word viewFactorWalls;


    // Private data

        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Map distributed
        autoPtr<IOdistributionMap> map_;

        //- Coarse mesh
        singleCellFvMesh coarseMesh_;

        //- Net radiative heat flux [W/m2]
        volScalarField qr_;

        //- View factor matrix
        autoPtr<scalarSquareMatrix> Fmatrix_;

        //- Inverse of C matrix
        autoPtr<scalarSquareMatrix> CLU_;

        //- Selected patches
        labelList selectedPatches_;

        //- Total global coarse faces
        label totalNCoarseFaces_;

        //- Total local coarse faces
        label nLocalCoarseFaces_;

        //- Constant emissivity
        bool constEmissivity_;

        //- Iterations Counter
        label iterCounter_;

        //- Pivot Indices for LU decomposition
        labelList pivotIndices_;

        //- Use Solar Load model
        bool applySolarLoad_;

        //- Solar load radiation model
        autoPtr<solarLoad> solarLoad_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Insert view factors into main matrix
        void insertMatrixElements
        (
            const globalIndex& index,
            const label fromProci,
            const labelListList& globalFaceFaces,
            const scalarListList& viewFactors,
            scalarSquareMatrix& matrix
        );


public:

    //- Runtime type information
    TypeName("viewFactor");


    // Constructors

        //- Construct from components
        viewFactor(const volScalarField& T, const dimensionedScalar& TRef);

        //- Construct from components
        viewFactor
        (
            const dictionary& dict,
            const volScalarField& T,
            const dimensionedScalar& TRef
        );

        //- Disallow default bitwise copy construction
        viewFactor(const viewFactor&) = delete;


    //- Destructor
    virtual ~viewFactor();


    // Member functions

        //- Solve system of equation(s)
        void calculate();

        //- Read radiation properties dictionary
        bool read();

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        virtual tmp<volScalarField::Internal> Ru() const;

        //- Const access to total radiative heat flux field
        inline const volScalarField& qr() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viewFactor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radiationModels/viewFactor/viewFactorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
