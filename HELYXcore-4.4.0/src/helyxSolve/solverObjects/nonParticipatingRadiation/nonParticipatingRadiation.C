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
    (c) 2010-2019 Engys Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "nonParticipatingRadiation.H"
#include "solverObjects/solverOption/SolverOption.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(nonParticipatingRadiation, 0);
}
}

makeFvSolverOption(nonParticipatingRadiation);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::nonParticipatingRadiation::nonParticipatingRadiation
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    Tname_(),
    radiation_()
{
    read(dict);

    Info<< endl;

    // must correct radiation before 1st temperature solve
    // otherwise Qin will be zero and the initial T boundary values will be too
    // low
    if (debug)
    {
        Info<< "    " << "Solving non-participating DOM radiation" <<  endl;
    }
    radiation_->correct();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::nonParticipatingRadiation::correct(const word& solveName)
{
    if (debug)
    {
        Info<< "    " << "Solving non-participating DOM radiation" <<  endl;
    }

    //check and modify if necessary fvSchemes and fvSolution settings
    radiation_->correct();
}

void Foam::fv::nonParticipatingRadiation::write()
{
    // implementation in wrapper class for functionObject-based solvers
    // no additional solver-specific write statement needed here
}


void Foam::fv::nonParticipatingRadiation::read(const dictionary& dict)
{
    Tname_ = dict.lookupOrDefault<word>("T", "T");
    dimensionedScalar TRef(basicThermo::TRefIfFound(obr_));
    radiation_ =
    (
        radiationModel::New
        (
            mesh_.lookupObject<volScalarField>(Tname_),
            TRef
        )
    );
}



// ************************************************************************* //
