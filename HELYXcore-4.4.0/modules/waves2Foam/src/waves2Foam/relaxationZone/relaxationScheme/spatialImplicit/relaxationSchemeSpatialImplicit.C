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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "relaxationSchemeSpatialImplicit.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationSchemeSpatialImplicit, 0);
addToRunTimeSelectionTable
(
    relaxationScheme,
    relaxationSchemeSpatialImplicit,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationSchemeSpatialImplicit::relaxationSchemeSpatialImplicit
(
    const word& subDictName,
    const fvMesh& mesh,
    vectorField& U,
    scalarField& alpha
)
:
    relaxationScheme(subDictName, mesh, U, alpha),
    factor_(coeffDict_.lookupOrDefault<scalar>("factor",1.0))
{
    const scalarField sigma(relaxShape_->sigma());

    weight_.setSize(sigma.size(), 0.0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


volScalarField& relaxationSchemeSpatialImplicit::relaxationWeightsMomentum()
{
    if (!mesh_.thisDb().foundObject<volScalarField>("relaxationWeightsMomentum"))
    {
        FatalErrorIn("relaxationSchemeSpatialImplicit::relaxationWeightsMomentum()")
            << "The field 'relaxationWeights' is not found. Use a solver, that"
            << "\nsupports the use of implicit relaxation."
            << endl << endl << exit(FatalError);
    }
    return const_cast<volScalarField&>
    (
        mesh_.thisDb().lookupObject<volScalarField>("relaxationWeightsMomentum")
    );
}


volScalarField& relaxationSchemeSpatialImplicit::targetAlpha()
{
    if (!mesh_.thisDb().foundObject<volScalarField>("targetAlphaField"))
    {
        FatalErrorIn("relaxationSchemeSpatialImplicit::targetAlpha()")
            << "The field 'targetAlphaField' is not found. Use a solver, which"
            << "\nsupport the use of implicit relaxation."
            << endl << endl << exit(FatalError);
    }
    return const_cast<volScalarField&>
    (
        mesh_.thisDb().lookupObject<volScalarField>("targetAlphaField")
    );
}


volVectorField& relaxationSchemeSpatialImplicit::targetVelocity()
{
    if (!mesh_.thisDb().foundObject<volVectorField>("targetVelocityField"))
    {
        FatalErrorIn("relaxationSchemeSpatialImplicit::targetVelocity()")
            << "The field 'targetVelocityField' is not found. Use a solver, that"
            << "\nsupports the use of implicit relaxation."
            << endl << endl << exit(FatalError);
    }
    return const_cast<volVectorField&>
    (
        mesh_.thisDb().lookupObject<volVectorField>("targetVelocityField")
    );
}


void relaxationSchemeSpatialImplicit::correct()
{
    // Obtain relaxation zone cells and local sigma coordinate
    // The number of cells and the sigma coordinate can have changed
    // for dynamic meshes
    const labelList& cells = relaxShape_->cells();
    const scalarField& sigma = relaxShape_->sigma();


    // Compute the relaxation weights - only changes for moving/changing meshes
    if (weight_.size() != sigma.size())
    {
        weight_.setSize( sigma.size(), 0.0 );
    }

    relaxWeight_->weights(cells, sigma, weight_);

    // Perform the correction
    const scalarField& V ( mesh_.V() );
    const vectorField& C ( mesh_.C() );
    const cellList& cc( mesh_.cells() );
    const pointField& pp( mesh_.points() );
    const faceList& fL( mesh_.faces() );
    // Get reference to implicit fields
    volScalarField& relaxWeights = relaxationWeightsMomentum();
    volScalarField& targetAlphaField = targetAlpha();
    volVectorField& targetVelocityField = targetVelocity();
    //relaxWeights *= 0.0;
    // Loop over the cells to relax
    forAll(cells, celli)
    {
        const label cellNo = cells[celli];
        const pointField& p = cc[cellNo].points(fL, pp);

        // Evaluate the cell height and the signedDistance to the surface from
        // the cell centre
        scalar cellHeight
            (
                Foam::max( p & waveProps_->returnDir() )
              - Foam::min( p & waveProps_->returnDir() )
            );
        scalar sdc( signedPointToSurfaceDistance( C[cellNo] ) );

        // Target variables
        scalar alphaTarget( 0.0 );
        vector UTarget( waveProps_->windVelocity( mesh_.time().value() ) );


        // Only do cutting, if surface is close by
        if (Foam::mag(sdc) <= 2*cellHeight)
        {

            localCellNeg lc(mesh_, cellNo);
            dividePolyhedral( point::zero, vector::one, lc);
            if (lc.ccNeg().size() >= 4)
            {
                UTarget = waveProps_->U(lc.centreNeg(), mesh_.time().value());
                alphaTarget = lc.magNeg()/V[cellNo];
            }

        }
        else if (sdc < 0)
        {
            alphaTarget = 1.0;
            UTarget     = waveProps_->U( C[cellNo], mesh_.time().value() );
        }
        else
        {
            alphaTarget = 0.0;
            UTarget     = waveProps_->windVelocity( mesh_.time().value() );
        }

        relaxWeights[cells[celli]] = factor_*(1.0 - weight_[celli]);
        targetAlphaField[cells[celli]] = alphaTarget;
        targetVelocityField[cells[celli]] = UTarget;

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// ************************************************************************* //
