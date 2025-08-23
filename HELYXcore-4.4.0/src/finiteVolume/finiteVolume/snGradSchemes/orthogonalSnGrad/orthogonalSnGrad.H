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
    (c) 2024 Engys Ltd.

Class
    Foam::fv::orthogonalSnGrad

Group
    grpFvSnGradSchemes

Description
    Simple central-difference snGrad scheme using the cell-centre to cell-centre
    delta-coefficients.

SourceFiles
    orthogonalSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef orthogonalSnGrad_H
#define orthogonalSnGrad_H

#include "finiteVolume/snGradSchemes/snGradScheme/snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class orthogonalSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class orthogonalSnGrad
:
    public snGradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("orthogonal");


    // Constructors

        //- Construct from mesh
        orthogonalSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}


        //- Construct from mesh and data stream
        orthogonalSnGrad
        (
            const fvMesh& mesh, Istream&, const word& forceGradSchemeName
        )
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~orthogonalSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the orthogonalSnGrad
        //  for the given field
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const orthogonalSnGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/snGradSchemes/orthogonalSnGrad/orthogonalSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
