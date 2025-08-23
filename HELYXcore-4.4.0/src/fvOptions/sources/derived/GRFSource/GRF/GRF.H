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
    (c) 2018-2024 Engys Ltd.

Class
    Foam::GRF

Description
    General reference frame class. This class keeps all the information
    about the axis, the origin, the rotational speed and faces that are part
    of this frame.


SourceFiles
    GRF.C

\*---------------------------------------------------------------------------*/


#ifndef GRF_H
#define GRF_H

#include "cellSetOption/cellSetOption.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "referenceFrames/coordinateFrame.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"
#include "sources/derived/frameSources/frameSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class GRF Declaration
\*---------------------------------------------------------------------------*/

class GRF
:
    public frameSources
{
protected:

    // Protected Data

        //- GRFDict
        dictionary GRFDict_;


    // Protected Member Functions

        //- Initialise
        void init();


public:

    // Constructors

        //- Construct from componets
        GRF
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        GRF(const GRF&) = delete;


    //- Destructor
    virtual ~GRF()
    {}

        //- Return equations for which we provide a source
        virtual void sourceFields(wordList& fieldNames);

        // Frame manipulations

        scalar getTriangularFlux(const face& f) const;

        void applyRotationalFlux(surfaceScalarField& phiOmega) const;

        void relaxBoundary
        (
            fvBlockMatrix<vector>& UEqn,
            const scalar
        ) const;

        void makeRelative(volVectorField& U) const;

        //- Filter-out the GRF region contribution from the given field
        // setting the corresponding values to zero
        virtual void zero(surfaceScalarField& phi) const;


        // I-O

            //- Read dictionary
            virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GRF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
