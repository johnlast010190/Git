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
    Foam::waveMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

SourceFiles
    waveMethod.C

\*---------------------------------------------------------------------------*/

#ifndef waveMethod_H
#define waveMethod_H

#include "meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class waveMethod Declaration
\*---------------------------------------------------------------------------*/

class waveMethod
:
    public meshToMeshMethod
{
public:

    //- Run-time type information
    TypeName("wave");


    // Constructors

        //- Construct from source and target meshes
        waveMethod
        (
            const polyMesh& src,
            const polyMesh& tgt,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        waveMethod(const waveMethod&) = delete;


    //- Destructor
    virtual ~waveMethod();


    // Member Functions

        //- Calculate addressing
        static void calculate
        (
            const polyMesh& src,
            const polyMesh& tgt,
            labelList& srcToTgtAddr
        );

        //- Calculate addressing and weights and optionally offset vectors
        virtual void calculate
        (
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght,
            pointListList& srcToTgtVec,
            labelListList& tgtToSrcAddr,
            scalarListList& tgtToSrcWght,
            pointListList& tgtToSrcVec
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const waveMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
