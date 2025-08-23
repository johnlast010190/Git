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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/error/error.H"

#include "processes/UOprocess/UOprocess.H"
#include "Kmesh/Kmesh.H"
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

complexVector UOprocess::WeinerProcess()
{
    return RootDeltaT*complexVector
    (
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>()),
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>()),
        complex(GaussGen.GaussNormal<scalar>(), GaussGen.GaussNormal<scalar>())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
UOprocess::UOprocess
(
    const Kmesh& kmesh,
    const scalar deltaT,
    const dictionary& UOdict
)
:
    GaussGen(),
    Mesh(kmesh),
    DeltaT(deltaT),
    RootDeltaT(sqrt(DeltaT)),
    UOfield(Mesh.size()),

    Alpha(UOdict.lookup<scalar>("UOalpha")),
    Sigma(UOdict.lookup<scalar>("UOsigma")),
    Kupper(UOdict.lookup<scalar>("UOKupper")),
    Klower(UOdict.lookup<scalar>("UOKlower")),
    Scale((Kupper - Klower)*pow(scalar(Mesh.size()), 1.0/vector::dim))
{
    const vectorField& K = Mesh;

    scalar sqrKupper = sqr(Kupper);
    scalar sqrKlower = sqr(Klower) + SMALL;
    scalar sqrK;

    forAll(UOfield, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrKupper && sqrK > sqrKlower)
        {
            UOfield[i] = Scale*Sigma*WeinerProcess();
        }
        else
        {
            UOfield[i] = complexVector
            (
                complex(0, 0),
                complex(0, 0),
                complex(0, 0)
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const complexVectorField& UOprocess::newField()
{
    const vectorField& K = Mesh;

    label count = 0;
    scalar sqrKupper = sqr(Kupper);
    scalar sqrKlower = sqr(Klower) + SMALL;
    scalar sqrK;

    forAll(UOfield, i)
    {
        if ((sqrK = magSqr(K[i])) < sqrKupper && sqrK > sqrKlower)
        {
            count++;
            UOfield[i] =
                (1.0 - Alpha*DeltaT)*UOfield[i]
              + Scale*Sigma*WeinerProcess();
        }
    }

    Info<< "    Number of forced K = " << count << nl;

    return UOfield;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
