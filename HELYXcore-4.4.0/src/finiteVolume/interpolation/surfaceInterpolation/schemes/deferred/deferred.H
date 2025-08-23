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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd

Class
    Foam::deferred

Group
    grpFvdeferreds

Description
    Deferred-correction scheme - turns any interpolation scheme into an upwind-
    weighted scheme with and explicit deferred correction.

SourceFiles
    deferred.C

\*---------------------------------------------------------------------------*/

#ifndef deferred_H
#define deferred_H

#include "interpolation/surfaceInterpolation/surfaceInterpolationScheme/surfaceInterpolationScheme.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class deferred Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class deferred
:
    public upwind<Type>
{
protected:

    // Protected data

        tmp<surfaceInterpolationScheme<Type>> scheme_;


public:

    //- Runtime type information
    TypeName("deferred");


    // Constructors

        //- Construct from mesh, faceFlux and scheme
        deferred
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const autoPtr<surfaceInterpolationScheme<Type>>& scheme
        )
        :
            upwind<Type>(mesh, faceFlux),
            scheme_(scheme.ptr())
        {}


        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        deferred
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& is
        )
        :
            upwind<Type>(mesh, db, is),
            scheme_(surfaceInterpolationScheme<Type>::New(mesh, db, is))
        {}


        //- Construct from mesh, flux and Istream.
        deferred
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux, is),
            scheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        deferred(const deferred&) = delete;


    //- Destructor
    virtual ~deferred() {};


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const override
        {
            return true;
        }

        //- Return explicit correction
        tmp<SurfaceField<Type>> correction
        (
            const VolField<Type>& vf
        ) const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const deferred&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/schemes/deferred/deferredTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
