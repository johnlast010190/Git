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
    (c) 2016 OpenCFD ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::maxDeltaxyz

Description
    Delta calculated by taking the maximum distance between the cell centre
    and any face centre.  For a regular hex cell, the computed delta will
    equate to half of the cell width; accordingly, the deltaCoeff model
    coefficient should be set to 2 for this case.

SourceFiles
    maxDeltaxyz.C

\*---------------------------------------------------------------------------*/

#ifndef maxDeltaxyzDelta_H
#define maxDeltaxyzDelta_H

#include "LES/LESdeltas/LESdelta/LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class maxDeltaxyz Declaration
\*---------------------------------------------------------------------------*/

class maxDeltaxyz
:
    public LESdelta
{
    // Private data

        //- Model coefficient
        scalar deltaCoeff_;


public:

    //- Runtime type information
    TypeName("maxDeltaxyz");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        maxDeltaxyz
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        maxDeltaxyz(const maxDeltaxyz&) = delete;


    //- Destructor
    virtual ~maxDeltaxyz()
    {}


    // Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const maxDeltaxyz&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
