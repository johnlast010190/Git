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
    (c) 2011-2024 Engys Ltd.

Class
    Foam::functionObjects::fieldBlendingFactor

Group
    grpFieldFunctionObjects

Description
    Calculates the field blending factors for use with localBlended scheme.
    A blending factor of 1 corresponds to pure scheme 1, while a blending factor
    of 0 corresponds to pure scheme 2 with in-between values constituting a
    linear combination. Least stable scheme should always be scheme 1.

Usage
    Example of function object specification:
    \verbatim
    fieldBlendingFactor1
    {
        // Type of functionObject
        type fieldBlendingFactor;

        // Where to load it from (if not already in solver)
        functionObjectLibs ("libfieldFunctionObjects.so");


        blendingSources // static blending sources
        (
            nonortho1
            {
                type nonOrthogonality;
                minValue 70; //non-orthogonality in degrees
                maxValue 85;
                minCoeff 0-1;
                maxCoeff 0-1;
                functionType linear|cubic;
            }
            faceSkewness
            maxFaceVolRatio
            maxCellVolume
            minCellVolume
            maxPointVolRatio
            cs1
            {
                type cellSet;
                posCoeff 1;
                interceptValue 0.5;
                negCoeff 0;
                functionType heaviside;
                //pos values use posCoeff, neg value use negCoeff
                //==interceptValue returns 0.5*(posCoeff+negCoeff)
                //cells that are part of the set have positive values
                //the rest have negative values
                setSources
                (
                    boxToCell
                    {
                       box   (0 0 0) (1 1 1);
                    }
                );
            }
            wd1
            {
                type wallDistance;
                minValue 70; //non-rthogonality in degrees
                maxValue 85;
                minCoeff 0-1;
                maxCoeff 0-1;
                functionType linear|cubic;
            }
            fs1
            {
                type faceSet;
                interceptValue 0.5;
                posCoeff 1;
                negCoeff 0;
                functionType heaviside;
                setSources
                (
                    boxToFace
                    {
                       box   (0 0 0) (1 1 1);
                    }
                );
            }


            des1
            {
                type DES;
                minValue 0.9; // (kappa*y) / (deltaCoeff*delta)
                maxValue 1.1;
                minCoeff 0;
                maxCoeff 1;
                functionType cubic;
                kappa 0.41;
                delta cubeRootVol; //LES cell length scale
                cubeRootVolCoeffs
                {
                    deltaCoeff 0.158;
                }
            }

        );

        fieldBlending
        (
            U
            {
                // limits relative increase in blending factor
                stabilise 10;

                blendingSources
                (
                    cfl1
                    {
                        type cfl;
                        minValue 2;
                        maxValue 5;
                        minCoeff 1;
                        maxCoeff 0;
                        functionType cubic;
                        phi phi;
                        rho rho; //compressible only
                        write   true;
                    }

                    Peclet1
                    {
                        type Pe;
                        minValue 10;
                        maxValue 20;
                        minCoeff 1;
                        maxCoeff 0.9;
                        functionType cubic;
                        phi phi;
                    }
                );
            }
        );
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldBlendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldBlendingFactor_H
#define functionObjects_fieldBlendingFactor_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fieldBlendingFactor/fieldBlendingItem/fieldBlendingItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldBlendingFactor Declaration
\*---------------------------------------------------------------------------*/

class fieldBlendingFactor
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- time index record
        label lastUpdateIndex_;

        // static blending mask
        autoPtr<fieldBlendingItem> staticBlend_;

        // PtrList of dynamic blending factor objects
        PtrList<fieldBlendingItem> fieldBlendingItems_;


public:

    //- Runtime type information
    TypeName("fieldBlendingFactor");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        fieldBlendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fieldBlendingFactor(const fieldBlendingFactor&) = delete;


    //- Destructor
    virtual ~fieldBlendingFactor();


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary&) override;

        //- Execute the blending factor calculation
        virtual bool execute() override;

        //- Calculate the field blending factors and write
        virtual bool write() override;

        bool end() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldBlendingFactor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
