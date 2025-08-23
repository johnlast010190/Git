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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2017-2024 Engys Ltd.

Class
    Foam::fv::anisotropicConductivityTransformation

Description
    Local coordinate transformation

    Sources described by:



        Anialpha
        {
            type            anisotropicConductivityTransformation;
            active          true;
            selectionMode   all;

            anisotropicConductivityTransformationCoeffs
            {

                coordinateSystem
                {
                    type cylindrical;

                    origin (0 0.1 0);
                    coordinateRotation
                    {
                        type    axesRotation;
                        e1  (1 0 0);
                        e3  (0 0 1);
                    }

                }
            }
        }

Note:
    The 'uniform' option allows the use of a time-varying uniform temperature
    by means of the DataEntry type.

SourceFiles
    fvOption.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicConductivityTransformation_H
#define anisotropicConductivityTransformation_H

#include "cellSetOption/cellSetOption.H"
#include "coordinate/systems/coordinateSystem.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class anisotropicConductivityTransformation Declaration
\*---------------------------------------------------------------------------*/

class anisotropicConductivityTransformation
:
    public cellSetOption
{

    // Private

        //- coordinate system
        autoPtr<coordinateSystem> csysPtr_;

        //- Reference Frame Ptr
        coordinateFrame* coorFramePtr_;

        //- local coordinate field
        autoPtr<volTensorField> localCartesian_;


    // Private Member Functions

        //- Local coordinate system
        const coordinateSystem& csys() const;

        //- Transform principal values of symmTensor
        symmTensor transformPrincipal
        (
            const tensor& tt,
            const vector& st
         ) const;

        //- Transform principal values of symmTensor
        void transformField
        (
            symmTensorField& fld,
            const tensorField& tt,
            const vectorField& st
        ) const;

        //- Transform principal values of symmTensor
        void transformField
        (
            volSymmTensorField& fld,
            const volTensorField& tt,
            const volVectorField& st
        ) const;

        //- Calculate local Cartesian coordinate
        tmp<volTensorField> updateTransformationTensor() const;


public:

    //- Runtime type information
    TypeName("anisotropicConductivityTransformation");


    // Constructors

        //- Construct from components
        anisotropicConductivityTransformation
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        anisotropicConductivityTransformation
        (
            const anisotropicConductivityTransformation&
        ) = delete;


    //- Destructor
    virtual ~anisotropicConductivityTransformation()
    {}


    // Member Functions

        virtual void sourceFields(wordList& fieldNames);

        virtual bool alwaysApply() const;


        // Correct

            //- Apply tensorial condutivity
            virtual void correct(volSymmTensorField& kappa);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const anisotropicConductivityTransformation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
