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
    (c) 2020 Engys Ltd.

Class
    Foam::roughWallCoefficients

Description
    Class that keeps information about the rough wall coefficients.

SourceFiles
    roughWallCoefficients.C

 \*---------------------------------------------------------------------------*/

#ifndef roughWallCoefficients_H
#define roughWallCoefficients_H

#include "fields/Fields/fieldTypes.H"
#include "fvMesh/fvPatches/fvPatch/fvPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class roughWallCoefficients Declaration
\*---------------------------------------------------------------------------*/
class roughWallCoefficients {

    scalar roughnessHeight_;

    scalar roughnessConstant_;

    scalar roughnessFactor_;

    void checkOldFormat(const dictionary& dict);

 public:

    //Constructors
    roughWallCoefficients();

    roughWallCoefficients(const dictionary &dict);


    bool roughnessIsActive() const;

    inline scalar roughnessHeight() const
    {
        return roughnessHeight_;
    }

    inline scalar roughnessConstant() const
    {
        return roughnessConstant_;
    }

    inline scalar roughnessFactor() const
    {
        return roughnessFactor_;
    }

    void writeLocalEntries(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //