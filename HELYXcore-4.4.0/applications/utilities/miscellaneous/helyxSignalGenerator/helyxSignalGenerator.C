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
    (c) 2025 Engys Ltd.

Application
    Create a random signal that follows a normal distribution with a given
    mean and standard deviation. A transient part can be added in the
    beginning. An auto-correlation coefficient can also be defined.


\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "cfdTools/general/include/fvCFD.H"

int main(int argc, char *argv[])
{
    Foam::Info<< "Create time\n" << Foam::endl;

    Foam::argList args(argc, argv);
    Foam::Time runTime(Foam::Time::controlDictName, args);

    Foam::IOdictionary signalDict
    (
        IOobject
        (
            "signalDict",
            runTime.system(),
            "",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    label samples = signalDict.lookupType<label>("nSamples");
    scalar mean = signalDict.lookupType<scalar>("mean");
    scalar stdDev = signalDict.lookupType<scalar>("stddev");
    label transientLength = signalDict.lookupType<label>("transientLength");
    scalar correlation = signalDict.lookupType<scalar>("correlation");

    Random rng(Foam::label(111));

    scalarField signal(samples);
    scalar previousValue = mean;
    for (int i=0; i<samples; ++i)
    {
            scalar noise = stdDev * rng.GaussNormal<scalar>();
            if (i == 0)
            {
                signal[i] = mean + noise;
            }
            else
            {
                signal[i] = correlation * previousValue + (1 - correlation) * (mean + noise);
//                signal[i] = correlation * previousValue + (mean + noise);
            }
            previousValue = signal[i];
    }

    scalar pi = 3.141592;
    if (transientLength>0)
    {
        scalar T = transientLength/2.;
        for (label i = 0; i < transientLength; ++i)
        {
//            signal[i] = mean + (signal[i] - mean) * scalar(i) / transientLength;
            signal[i] += 0.005*Foam::cbrt(static_cast<scalar>(i))*Foam::sin
            (static_cast<scalar>(2*pi*i/(T)));
        }
    }

    std::ofstream outFile("signal.dat");

    outFile <<"#Time"<<tab<<"drag"<<nl;
    for (label i = 0; i < samples; ++i)
    {
        outFile << i+1 << tab << signal[i] << "\n";
    }
    outFile.close();

    Foam::Info<< "Signal written to file " << Foam::endl;

    return 0;
}
