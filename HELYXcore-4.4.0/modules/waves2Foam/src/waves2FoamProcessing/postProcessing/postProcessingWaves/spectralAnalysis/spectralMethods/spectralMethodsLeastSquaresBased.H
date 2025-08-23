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
    Foam::spectralMethodsLeastSquaresBased

Description
    Helper functions for spectral analysis using a least-squares based
    approach.

SourceFiles
    spectralMethodsLeastSquaresBased.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef spectralMethodsLeastSquaresBased_H
#define spectralMethodsLeastSquaresBased_H

#include "cfdTools/general/include/fvCFD.H"
#include "complexExp.H"

#if EXTBRANCH==1
    #include "scalarSquareMatrix.H"
#elif OFPLUSBRANCH==1
    #include "matrices/scalarMatrices/scalarMatrices.H"
#else
    #include "matrices/scalarMatrices/scalarMatrices.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class spectralMethodsLeastSquaresBased
{
private:

    // Private Member Functions

        void computePowerSpectrum
        (
            const scalarField&,
            const scalarField&,
            const label&,
            const scalar&,
            scalarField&
        );


public:

    //- Runtime type information
    TypeName("spectralMethodsLeastSquaresBased");


    // Constructors

        //- Construct from components
        spectralMethodsLeastSquaresBased
        (
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        spectralMethodsLeastSquaresBased
        (
            const spectralMethodsLeastSquaresBased&
        ) = delete;


    //- Destructor
    virtual ~spectralMethodsLeastSquaresBased();


    // Member Functions

        //- Solves a least squares problem, with the columns given in the
        //  List<scalarField> and the right hand side in the scalarField.
        //  The solution is returned in the right hand side field
        void solve
        (
            const List<scalarField>&,
            scalarField&
        );

        scalarField frequencies
        (
            const label&
        );

        scalarField powerSpectra
        (
               const scalarField&,
            const scalarField&,
            const label&,
            const scalar&
        );

        List<vectorField> powerSpectra
        (
               const scalarField&,
            const List<vectorField>&,
            const label&,
            const scalar&
        );

        List<scalarField> powerSpectra
        (
               const scalarField&,
            const List<scalarField>&,
            const label& N,
            const scalar&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const spectralMethodsLeastSquaresBased&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
