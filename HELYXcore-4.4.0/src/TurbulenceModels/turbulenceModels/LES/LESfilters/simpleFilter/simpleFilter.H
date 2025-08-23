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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::simpleFilter

Description
    Simple top-hat filter used in dynamic LES models.

    Implemented as a surface integral of the face interpolate of the field.

SourceFiles
    simpleFilter.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFilter_H
#define simpleFilter_H

#include "LES/LESfilters/LESfilter/LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleFilter Declaration
\*---------------------------------------------------------------------------*/

class simpleFilter
:
    public LESfilter

{
public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct from components
        simpleFilter(const fvMesh& mesh);

        //- Construct from IOdictionary
        simpleFilter(const fvMesh& mesh, const dictionary&);

        //- Disallow default bitwise copy construction
        simpleFilter(const simpleFilter&) = delete;


    //- Destructor
    virtual ~simpleFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const simpleFilter&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
