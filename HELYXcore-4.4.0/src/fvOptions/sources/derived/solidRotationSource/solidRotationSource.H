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
    (c) 2012-2013 OpenFOAM Foundation

Class
    Foam::fv::solidRotationSource

Description

    Supplies rotational convective heat flux for solid equation, using GRF

    This source is primarily intended to be used in zones with rotational
    symmetry, allowing rotation of solids to be simulated without mesh movement.

    However, deviations from rotational symmetry are tolerated by assuming no
    flux into/out of the boundaries of the zone, and using the assumption of
    volumetric conservation to modify the rotational flux to suit this
    constraint. This means that the rotational convection will be blocked by
    any gaps in the mesh, but will flow around these were possible, effectively
    yielding a useful approximation for solids that are not bodies of
    revolution.

    This is achieved by solving for a variable "rotatingFluxPotential", which
    needs solver information in fvSolution and a
    laplacian(rotatingFluxPotential) scheme specified in fvSchemes. Defaults
    will automatically be applied if not specified.

    \heading Source usage

    Example usage:
    \verbatim
    disc_MRF
    {
        type            solidRotationSource;
        active          true;
        referenceFrames
        {
            all
            {
                cellZone        disk;
                selectionMode   cellZone;
                referenceFrame  discFrame;
            }
        }
        fields          ( T );

        // Number of non-orthogonal correctors to use in correcting the
        // rotational flux. Default = 0
        nNonOrthogonalCorrectors 1;

        // Whether to write out a field for post-processing called
        // 'rotationalVelocity' to visualise the corrected convection flow
        // that is calculated for non-rotationally symmetric bodies.
        // Default = false
        writeRotationalVelcity false;
    }
    \endverbatim

    Convection scheme may be specified as
    \verbatim
    divSchemes
    {
        div(phiOmega,h) Gauss linearUpwind limitedGrad(h);
    }
    \endverbatim
    In addition, a cell-limited gradient scheme is advisable
    \verbatim
    gradSchemes
    {
        limitedGrad(h) cellLimited Gauss linear 1;
    }
    \verbatim

    Do not use 'bounded' as this is accounted for explicitly.

    The above schemes are automatically added as defaults, so they need not
    be specified unless changing the schemes.

    See also: GRFSource

SourceFiles
    solidRotationSource.C

\*---------------------------------------------------------------------------*/
#ifndef solidRotationSource_H
#define solidRotationSource_H

#include "sources/derived/GRFSource/GRFSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class solidRotationSource Declaration
\*---------------------------------------------------------------------------*/

class solidRotationSource
:
    public GRFSource
{
protected:

    // Protected Member Variables

        //- Rotating flux correction potential
        mutable volScalarField p_;

        //- Rotational flux
        mutable autoPtr<surfaceScalarField> phi_;

        //- Rotational velocity, if requested for postprocessing
        mutable autoPtr<volVectorField> rotationalVelocity_;

        //- The uncorrected GRF flux used to produce the current phi
        //  Stored so it can be checked for changes
        mutable surfaceScalarField phi0_;

        //- Number of non-orthogonal correctors
        label nNonOrthCorr_;

        //- Write out the rotational velocity reconstructed from the corrected
        //  flux, for post-processing purposes
        bool writeRotationalVelocity_;


    // Protected Member Functions

        virtual void addLocalFvSchemes(const wordList& fieldNames) override;

public:

    //- Runtime type information
    TypeName("solidRotationSource");


    // Constructors

        //- Construct from components
        solidRotationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        solidRotationSource(const solidRotationSource&) = delete;


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict) override;

        //- Return rotation induced flux field for the entire domain
        virtual tmp<surfaceScalarField> phiOmega() const override;

        //- Calculate the volumetric flux by correcting the rotational flux
        virtual void calculatePhi() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidRotationSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
