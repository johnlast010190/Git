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

Class
    Foam::UpwindFitData

Description
    Data for the quadratic fit correction interpolation scheme to be used with
    upwind biased stencil.
    - linearCorrection = true : fit calculated for corrected linear scheme
    - linearCorrection = false : fit calculated for corrected upwind scheme

SourceFiles
    UpwindFitData.C

\*---------------------------------------------------------------------------*/

#ifndef UpwindFitData_H
#define UpwindFitData_H

#include "interpolation/surfaceInterpolation/schemes/FitData/FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedUpwindCellToFaceStencil;

/*---------------------------------------------------------------------------*\
                    Class UpwindFitData Declaration
\*---------------------------------------------------------------------------*/

template<class Polynomial>
class UpwindFitData
:
    public FitData
    <
        UpwindFitData<Polynomial>,
        extendedUpwindCellToFaceStencil,
        Polynomial
    >
{
    // Private data

        //- For each face of the mesh store the coefficients to multiply the
        //  stencil cell values by if the flow is from the owner
        List<scalarList> owncoeffs_;

        //- For each face of the mesh store the coefficients to multiply the
        //  stencil cell values by if the flow is from the neighbour
        List<scalarList> neicoeffs_;


    // Private Member Functions

        //- Calculate the fit for the all the mesh faces
        //  and set the coefficients
        void calcFit();


public:

    TypeName("UpwindFitData");


    // Constructors

        //- Construct from components
        UpwindFitData
        (
            const fvMesh& mesh,
            const extendedUpwindCellToFaceStencil& stencil,
            const bool linearCorrection,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~UpwindFitData()
    {}


    // Member functions

        //- Return reference to owner fit coefficients
        const List<scalarList>& owncoeffs() const
        {
            return owncoeffs_;
        }

        //- Return reference to neighbour fit coefficients
        const List<scalarList>& neicoeffs() const
        {
            return neicoeffs_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/surfaceInterpolation/schemes/UpwindFitScheme/UpwindFitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
