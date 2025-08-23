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
    Foam::relaxationZone

  Description
    Base class for relaxationSchemes.
    (This is the class formerly known as "waveRelaxationOnEcstacy".)

  SourceFiles
      relaxationZone.C

  Author
  Niels Gjoel Jacobsen and Bo Terp Paulsen, The Technical University of Denmark
  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef relaxationZone_H
#define relaxationZone_H

#include "relaxationZone/relaxationScheme/relaxationScheme.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class relaxationZone Declaration
\*---------------------------------------------------------------------------*/

class relaxationZone
{
private:

    const fvMesh& mesh_;
    volVectorField& U_;
    volScalarField& alpha_;

    wordList relaxNames_;

    List<autoPtr<relaxationSchemes::relaxationScheme>> relaxSchemePtr_;

    void resetTargetFields();

    void correctBoundaries();

protected:

    mutable volScalarField* relaxationWeightsMomentum_;

    mutable volScalarField* relaxationWeightsSource_;

    mutable volScalarField* targetAlpha_;

    mutable volVectorField* targetVelocity_;

public:
     //- Runtime type information
    TypeName("relaxationZone");

    // Constructors

        //- Construct from components
        relaxationZone
        (
            const fvMesh&,
            volVectorField&,
            volScalarField&
        );

    // Destructor

        virtual ~relaxationZone()
        {}

    // Member Functions

        virtual void correct();

        tmp<volScalarField> numericalBeach();

    // Access Functions
        const volScalarField& relaxationWeightsMomentum() const;

        const volScalarField& relaxationWeightsSource() const;

        const volScalarField& targetAlphaField() const;

        const volVectorField& targetVelocityField() const;

};


} // End namespace Foam

#endif
