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
    (c) 2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::Curle

Group
    grpFieldFunctionObjects

Description
    Calculates the acoustic pressure based on Curle's analogy.

    Curle's analogy is implemented as:

    \f[
        p' = \frac{1}{4 \pi c_0}\frac{\vec d}{|\vec d|^2}\frac{d(F)}{d(t)}
    \f]


    where
    \vartable
        p'          | Curle's acoustic pressure [Pa] or [Pa (m3/rho)]
        c_0         | Reference speed of sound [m/s]
        \vec d      | Distance vector to observer locations [m]
        F           | Force [N] or [N (m3/rho)]
    \endvartable

Note
    Only the normal-pressure force is included in the force calculation


Usage
    Example of function object specification:
    \verbatim
    Curle1
    {
        type        Curle;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (surface1 surface2);
        c0          330;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                  | Required | Default value
        type         | Type name: Curle             | yes      |
        field        | Pressure field name          | no       | p
        result       | Acoustic pressure field name | no       | Curle
        patches      | Sound generation patch names | yes      |
        c0           | Reference speed of sound     | yes      |
    \endtable


See also
  - Foam::functionObjects::fieldExpression
  - Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    Curle.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Curle_H
#define functionObjects_Curle_H

#include "fieldExpression/fieldExpression.H"
#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "dimensionedTypes/dimensionedVector/dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class Curle Declaration
\*---------------------------------------------------------------------------*/

class Curle
:
    public fieldExpression
{
    // Private data

        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Area-averaged centre of patch faces
            dimensionedVector x0_;

            //- Reference speed of souund
            dimensionedScalar c0_;


protected:

        //- Calculate the acoustic pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("Curle");


    // Constructors

        //- Construct from Time and dictionary
        Curle
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        Curle(const Curle&) = delete;


    //- Destructor
    virtual ~Curle();


    // Member Functions

        //- Read the Curle data
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Curle&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
