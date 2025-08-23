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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2021 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistrySolver/EulerImplicit/EulerImplicit.H"
#include "fields/Fields/Field/SubField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::EulerImplicit
(
    const fluidMulticomponentThermo& thermo
)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("EulerImplicitCoeffs")),
    cTauChem_(coeffsDict_.lookup<scalar>("cTauChem")),
    cTp_(this->nEqns()),
    R_(this->nEqns()),
    J_(this->nEqns()),
    E_(this->nEqns() - 2)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::EulerImplicit<ChemistryModel>::~EulerImplicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::EulerImplicit<ChemistryModel>::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    const label nSpecie = this->nSpecie();

    // Map the composition, temperature and pressure into cTp
    for (int i=0; i<nSpecie; i++)
    {
        cTp_[i] = max(0, c[i]);
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie + 1] = p;

    // Calculate the reaction rate and Jacobian
    this->jacobian(0, cTp_, li, R_, J_);

    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;
    const scalar cTot = sum(c);

    for (label i=0; i<nSpecie; i++)
    {
        if (R_[i] < -SMALL)
        {
            tMin = min(tMin, -(cTp_[i] + SMALL)/R_[i]);
        }
        else
        {
            tMin = min
            (
                tMin,
                max(cTot - cTp_[i], 1e-5)/max(R_[i], SMALL)
            );
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Assemble the Euler implicit matrix for the composition
    scalarField& source = E_.source();
    for (label i=0; i<nSpecie; i++)
    {
        E_(i, i) = 1/deltaT - J_(i, i);
        source[i] = R_[i] + E_(i, i)*cTp_[i];

        for (label j=0; j<nSpecie; j++)
        {
            if (i != j)
            {
                E_(i, j) = -J_(i, j);
                source[i] += E_(i, j)*cTp_[j];
            }
        }
    }

    // Solve for the new composition
    scalarField::subField(cTp_, nSpecie) = E_.LUsolve();

    // Limit the composition and transfer back into c
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0, cTp_[i]);
    }

    // Euler explicit integrate the temperature.
    // Separating the integration of temperature from composition
    // is significantly more stable for exothermic systems
    T += deltaT*R_[nSpecie];
}


// ************************************************************************* //
