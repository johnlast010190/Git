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
    (c) 2017-2024 Engys Ltd.

Class
    GIBAreaSmoothing

Description
    Surface area implicit smoothing for GIB.

SourceFiles
    GIBAreaSmoothing.C

\*---------------------------------------------------------------------------*/

#ifndef GIBAreaSmoothing_H
#define GIBAreaSmoothing_H

#include "fvMesh/fvMesh.H"
#include "faMesh/geodeticWallDist/geodeticWallDist/geodeticWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class GIBAreaSmoothing Declaration
\*---------------------------------------------------------------------------*/

class GIBAreaSmoothing
{
    // Private Data

        const fvMesh& mesh_;

        faMesh areaMesh_;

        const label& patchID_;

        boolList constraintFaces_;

        dimensionedScalar r_;

        scalar standardDeviation_;

        scalar cutOffAngle_;

        //- Create zero band for smooth transition
        bool smoothTransition_;

        const volScalarField& G_;

        mutable areaScalarField* Gsm_;

        mutable areaScalarField* Gsource_;

        mutable edgeScalarField* DsF_;

        mutable volScalarField* Gsmoothed_;

        mutable geodeticWallDist* y_;


    // Private Member Functions

        labelList labelListPatch(const label&);

        void smooth();

        void createFAMesh();

        void createFAFields();

        labelList findIsolatedAreas();

        void setSmoothingOptions(const dictionary&);

        //- Clear all demand-driven data
        virtual void clearOutData();


public:

    //- Runtime type information
    TypeName("GIBAreaSmoothing");


    // Constructors

        //- Construct from components
        GIBAreaSmoothing
        (
            const fvMesh&,
            const dictionary&,
            const label&,
            boolList
        );


    //- Destructor
    virtual ~GIBAreaSmoothing()
    {
        clearOutData();
    };


    // Member Functions

        //- Update
        void update();

        tmp<scalarField> smoothSens() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
