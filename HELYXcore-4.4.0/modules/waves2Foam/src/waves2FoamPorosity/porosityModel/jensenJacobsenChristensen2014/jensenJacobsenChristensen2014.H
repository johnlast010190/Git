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

Class
    Foam::jensenJacobsenChristensen2014

Description


SourceFiles
    jensenJacobsenChristensen2014.C

Author
    Niels G. Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef jensenJacobsenChristensen2014_H
#define jensenJacobsenChristensen2014_H

#include "porosityModel/wavesPorosityModel.H"
#include "jjc2014Zones.H"
//#include "porosityZones/porosityZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class jensenJacobsenChristensen2014 Declaration
\*---------------------------------------------------------------------------*/

class jensenJacobsenChristensen2014
:
    public wavesPorosityModel
{

protected:

    // Protected data
//        porosityZones pZones_;
        jjc2014Zones pZones_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("jensenJacobsenChristensen2014");

    // Constructors

        //- from components
        jensenJacobsenChristensen2014
        (
            const fvMesh& mesh
        );


    // Destructor

        ~jensenJacobsenChristensen2014()
        {}


    // Member Functions

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            VolField<vector>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const geometricOneField&,
            VolField<vector>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const dimensionedScalar&,
            VolField<vector>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector>> ddt
        (
            const volScalarField&,
            VolField<vector>&
        );

        //- Update the porosity field
        virtual void updatePorosity();

        //- Return a volScalarField with the porosities
//        virtual tmp<volScalarField> porosity() const;
        const volScalarField& porosity() const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        virtual void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
