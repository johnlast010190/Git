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

#include "stokesFifthProperties.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesFifthProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    stokesFifthProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


scalar stokesFifthProperties::eval(scalar& k)
{
    scalar S  = 1.0/Foam::cosh(2*k * depth_);
    scalar C0 = Foam::sqrt(Foam::tanh(k*depth_));
    scalar C2 = Foam::sqrt(Foam::tanh(k*depth_))*(2.0 + 7.0*Foam::pow(S, 2.0))
        /( 4.0*Foam::pow(1.0 - S, 2.0) );
    scalar C4 = Foam::sqrt(Foam::tanh(k*depth_))
        *(4.0 + 32.0*S - 116.0*Foam::pow(S, 2.0) - 400.0*Foam::pow(S, 3.0)
        - 71.0*Foam::pow(S, 4.0) + 146.0*Foam::pow(S, 5.0))
        /(32.0*Foam::pow(1.0 - S, 5.0) );
    scalar D2 = - Foam::sqrt(1.0/Foam::tanh(k*depth_))/2.0;
    scalar D4 = Foam::sqrt(1.0/Foam::tanh(k*depth_))
        *(2.0 + 4.0*S + Foam::pow(S,2.0) + 2.0*Foam::pow(S,3.0))
        /(8.0*Foam::pow(1.0 - S,3.0));

    return Foam::sqrt(k/G_)*Q_ - 2.0*PI_/(period_*Foam::sqrt(G_*k))
         + C0 + Foam::pow((k*height_/2.0),2.0)*(C2 + D2/(k*depth_))
         + Foam::pow(k*height_/2.0, 4.0)*(C4 + D4/(k*depth_));
}


scalar stokesFifthProperties::waveNumber()
{
    scalar lower(1.0e-7);

    scalar upper = Foam::max
        (
            4.0*PI_/( period_*Foam::sqrt( Foam::mag(G_)*depth_)),
            2.0*PI_/( Foam::pow( period_, 2.0))
        );

    scalar middle(0.5*(lower + upper) );

    scalar valLower( eval( lower ) ),
           //valUpper( eval( upper ) ),
           valMiddle( eval( middle ) );

    while (true)
    {
        if (Foam::sign(valLower) == Foam::sign( valMiddle ))
        {
            lower    = middle;
            valLower = valMiddle;
        }
        else
        {
            upper    = middle;
            //valUpper = valMiddle;
        }

        middle = 0.5*( lower + upper );

        valMiddle = eval(middle);

        if (Foam::mag(valMiddle) < 1.0e-13)
        {
            break;
        }
    }

    return middle;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesFifthProperties::stokesFifthProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write)
{
    Info<< "\nConstructing: " << this->type() << endl;

    period_ = dict.lookup<scalar>("period");
    depth_  = dict.lookup<scalar>("depth");
    height_ = dict.lookup<scalar>("height");
    Q_      = dict.lookup<scalar>("stokesDrift");

    omega_  = 2.0*PI_/period_ ;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void stokesFifthProperties::set( Ostream& os)
{
    scalar k = waveNumber();

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    writeGiven( os, "height");
    writeGiven( os, "period" );

    writeGiven( os, "depth" );
    writeGiven( os, "stokesDrift");
    writeGiven( os, "direction" );

    if (dict_.found("Tsoft"))
    {
        writeGiven( os, "Tsoft");
    }

    writeGiven( os, "phi");

    if (write_)
    {
        vector direction(dict_.lookup<vector>("direction"));
        direction /= Foam::mag(direction);
        direction *= k;

        writeDerived(os, "waveNumber", direction);
        writeDerived(os, "omega", omega_);
    }

    // Write the relaxation zone
    writeRelaxationZone( os );

    // Write the closing bracket
    writeEnding( os );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
