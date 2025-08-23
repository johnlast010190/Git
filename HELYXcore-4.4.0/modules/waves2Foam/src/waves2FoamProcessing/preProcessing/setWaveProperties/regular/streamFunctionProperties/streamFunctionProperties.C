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

#include "streamFunctionProperties.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(streamFunctionProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    streamFunctionProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


streamFunctionProperties::streamFunctionProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write),

    localRT_(rT)
{
    Info<< "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void streamFunctionProperties::writeInputFile() const
{
    // Get the input information and write input file
    autoPtr<OFstream> inputFile;
    inputFile.reset(new OFstream("fenton.inp"));

    inputFile().precision(14);

    // Wave height and depth information
    scalar height = dict_.lookup<scalar>("height");
    scalar depth = dict_.lookup<scalar>("depth");

    inputFile() << "finite " << height/depth << endl;

    // Specify with period scale: period/wave length
    Switch specifyPeriod(dict_.lookup("specifyPeriod"));
    scalar    periodScale = 0;
    if (specifyPeriod)
    {
        periodScale = dict_.lookup<scalar>("period");

        inputFile() << "period " << height/(G_*Foam::sqr(periodScale)) << endl;
    }
    else
    {
        periodScale = dict_.lookup<scalar>("waveLength");

        inputFile() << "wavelength " << height/periodScale << endl;
    }

    // Mean flow scale
    Switch specifyEuler(dict_.lookup("specifyEuler"));
    scalar currentScale = 0;
    if (specifyEuler)
    {
        currentScale = dict_.lookup<scalar>("eulerVelocity");

        inputFile() << "Euler " << currentScale/Foam::sqrt(height*G_) << endl;
    }
    else
    {
        currentScale = dict_.lookup<scalar>("stokesVelocity");

        inputFile() << "Stokes " << currentScale/Foam::sqrt(height*G_) << endl;
    }

    // Get the controls
    label nModes(dict_.lookup<label>("N"));
    label nIter(dict_.lookup<label>("Niter"));

    inputFile() << nModes << " " << nIter;

    // Done writing the input file
}


void streamFunctionProperties::set( Ostream& os)
{
    // Write the input file to the ThirdParty fenton program
    writeInputFile();

    // Execute the fortran program 'fenton4Foam'. Name change not to collide
    // with other installations of the same code.
    system("fenton4Foam > /dev/null");

    // Read the output file
    IOdictionary fentonCoeffs
    (
        IOobject
        (
            "fenton4Foam",
            localRT_.constant(),
            localRT_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read the various output variables
    dimensionedScalar depthNonDim(fentonCoeffs.lookup("depth"));
    dimensionedScalar periodNonDim(fentonCoeffs.lookup("period"));
    dimensionedScalar uBarNonDim(fentonCoeffs.lookup("uBar"));
    label N = fentonCoeffs.lookup<label>("nModes");

    scalarField coeffsA("aCoeffs", fentonCoeffs, N);
    scalarField coeffsB("bCoeffs", fentonCoeffs, N);

    // Convert to the correct OpenFoam (dimensioned) format
    scalar depth = dict_.lookup<scalar>("depth");
    scalar waveNumber = depthNonDim.value()/depth;

    scalar period = periodNonDim.value()/Foam::sqrt(G_*waveNumber);

    scalar uBar = uBarNonDim.value()*Foam::sqrt(G_/waveNumber);
    coeffsA /= waveNumber;

    forAll(coeffsB, coeffi)
    {
        coeffsB[coeffi] *= (coeffi + 1)*Foam::sqrt(G_/waveNumber);
    }

    label oldPrecision = os.precision(14);

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven(os, "waveType");
    writeGiven(os, "N");

    if (dict_.found("Tsoft"))
    {
        writeGiven( os, "Tsoft");
    }

    writeGiven(os, "depth");
    writeGiven(os, "phi");
    writeGiven(os, "height");
    os() << endl;
    writeGiven(os, "specifyPeriod");
    writeDerived(os, "period", period);

    if (!Switch(dict_.lookup("specifyPeriod")))
    {
        writeGiven(os, "waveLength");
    }
    os() << endl;

    writeGiven(os, "specifyEuler");

    if (Switch(dict_.lookup("specifyEuler")))
    {
        writeGiven(os, "eulerVelocity");
    }
    else
    {
         writeGiven(os, "stokesVelocity");
    }
    os() << endl;

    // This part should compute the properties for stream function wave theory
    vector direction( vector(dict_.lookup("direction")));
    direction /= Foam::mag(direction);

    writeDerived(os, "omega", 2*M_PI/period);
    writeDerived(os, "waveNumber", waveNumber*direction);
    writeDerived(os, "uBar", uBar);
    writeDerived(os, "A", coeffsA);
    writeDerived(os, "B", coeffsB);

    os.precision(oldPrecision);

    // Write the relaxation zone
    writeRelaxationZone(os);

    // Write the closing bracket
    writeEnding( os );

    // Clean-up after the execution of the external programme
    if (!dict_.lookupOrDefault<Switch>("keepTemporaryFiles", false))
    {
        Foam::rm(localRT_.constant() + "/../fenton.inp");
        Foam::rm(localRT_.constant() + "/../fenton.out");
        Foam::rm(localRT_.constant() + "/fenton4Foam");
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
