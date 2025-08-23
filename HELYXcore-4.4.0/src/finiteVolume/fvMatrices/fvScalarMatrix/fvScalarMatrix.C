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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMatrices/fvScalarMatrix/fvScalarMatrix.H"
#include "fields/fvPatchFields/basic/extrapolatedCalculated/extrapolatedCalculatedFvPatchFields.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvMatrix<Foam::scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
               *value;
        }
    }
}


template<>
Foam::autoPtr<Foam::fvMatrix<Foam::scalar>::fvSolver>
Foam::fvMatrix<Foam::scalar>::solver
(
    const dictionary& solverControls
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    word regionName;
    if (psi_.mesh().name() != polyMesh::defaultRegion)
    {
        regionName = psi_.mesh().name() + "::";
    }
    addProfiling(solve, "fvMatrix::solve." + regionName + psi_.name());
    #endif

    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solver(const dictionary& solverControls) : "
               "solver for fvMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    autoPtr<fvMatrix<scalar>::fvSolver> solverPtr
    (
        new fvMatrix<scalar>::fvSolver
        (
            *this,
            lduMatrix::solver::New
            (
                psi_.name(),
                *this,
                boundaryCoeffs_,
                internalCoeffs_,
                psi_.boundaryField().scalarInterfaces(),
                solverControls
            )
        )
    );

    diag() = saveDiag;

    return solverPtr;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::fvSolver::solve
(
    const dictionary& solverControls
)
{
    VolField<scalar>& psi = const_cast<VolField<scalar>&>(fvMat_.psi());

    scalarField saveDiag(fvMat_.diag());
    fvMat_.addBoundaryDiag(fvMat_.diag(), 0);

    scalarField totalSource(fvMat_.source());
    fvMat_.addBoundarySource(totalSource, fvMat_.psi(), false);

    // Assign new solver controls
    solver_->read(solverControls);

    solverPerformance solverPerf = solver_->solve
    (
        psi.primitiveFieldRef(),
        totalSource
    );

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(fvMat_.mesh().comm()));
    }

    fvMat_.diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    VolField<scalar>& psi = const_cast<VolField<scalar>&>(psi_);

    this->boundaryManipulate(psi.boundaryFieldRef());

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    scalarField totalSource(source_);
    addBoundarySource(totalSource, psi_, false);

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        psi_.boundaryField().scalarInterfaces(),
        solverControls
    )->solve(psi.primitiveFieldRef(), totalSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Foam::scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.primitiveField(),
            source_ - boundaryDiag*psi_.primitiveField(),
            boundaryCoeffs_,
            psi_.boundaryField().scalarInterfaces(),
            0
        )
    );

    addBoundarySource(tres.ref(), psi_);

    return tres;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H() const
{
    tmp<volScalarField> tHphi
    (
        volScalarField::New
        (
            "H(" + psi_.name() + ')',
            psi_.db(),
            psi_.mesh(),
            dimensions_/dimVol,
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& Hphi = tHphi.ref();

    Hphi.primitiveFieldRef() = (lduMatrix::H(psi_.primitiveField()) + source_);
    addBoundarySource(Hphi.primitiveFieldRef(), psi_);

    Hphi.primitiveFieldRef() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H1() const
{
    tmp<volScalarField> tH1
    (
        volScalarField::New
        (
            "H(1)",
            psi_.db(),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1.ref();

    H1_.primitiveFieldRef() = lduMatrix::H1();
    //addBoundarySource(Hphi.primitiveField(), psi_);

    H1_.primitiveFieldRef() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


// ************************************************************************* //
