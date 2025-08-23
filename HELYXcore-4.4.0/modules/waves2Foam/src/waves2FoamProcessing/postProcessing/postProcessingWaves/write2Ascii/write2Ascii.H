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
    (c) 2024 Engys Ltd.

Class
    Foam::write2Ascii

Description
    Reads the OF-formatted files and write the content in ASCII format without
    any headers, etc.

SourceFiles
    write2Ascii.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef write2Ascii_H
#define write2Ascii_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class write2Ascii
:
    public postProcessingWaves
{
protected:

    // Protected data

        word inputDir_;


public:

    //- Runtime type information
    TypeName("write2Ascii");


    // Constructors

        //- Construct from components
        write2Ascii
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        write2Ascii(const write2Ascii&) = delete;


    //- Destructor
    virtual ~write2Ascii();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const write2Ascii&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
