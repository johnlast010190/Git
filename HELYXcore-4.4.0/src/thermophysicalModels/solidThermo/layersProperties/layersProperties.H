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
    (c) 2019 Engys Ltd.

Class
    Foam::layersProperties

Description
    Container class for solid layers material properties

    @verbatim
        wall
        {
            layers
            {
              layer1
              {
                  thickness 0.003;
                  kappa     16;
                  rcontact  0.0001;
              }
              layer2
              {
                  thickness 0.05;
                  kappa     10;
                  rcontact  0.00005;
              }
              layer3
              {
                  thickness 0.002;
                  kappa     100;
                  rcontact  0.0002;
              }
            }
        }
    @endverbatim

SourceFiles
    layersProperties.C

\*---------------------------------------------------------------------------*/

#ifndef layersProperties_H
#define layersProperties_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class layersProperties Declaration
\*---------------------------------------------------------------------------*/

class layersProperties
{

//private:

    //- layer name
    List<word> layerNames_;

    //- layer thickness [m]
    scalarList t_;

    //- layer thermal conductivity [W/m K]
    scalarList kappa_;

    //- layer contact resistance
    //  per unit of area [m^2 K/ W]
    scalarList rcontact_;

    //- total layers resistance
    //  per unit of area [m^2 K/ W]
    mutable scalar rtot_;

    // Private Member Functions
    void calculateTotalResistance() const;

public:

    // Dictionary constructor
    layersProperties
    (
        const dictionary& dict
    );

    // Construct as copy
    layersProperties
    (
        const layersProperties&
    );

    // Null constructor
    layersProperties(){}

    // Destructor
    virtual ~layersProperties(){}

    // Return total resistance per unit of area
    inline scalar totalResistance() const
    {
        return rtot_;
    }

    //- Read patch dictionary
    virtual bool read(const dictionary& dict);

    //- Write function
    virtual void write(Ostream& os) const;

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
