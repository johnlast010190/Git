/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2024 Engys Ltd.

Class
    solvePDECommonFiniteVolume

Description
    Common functionality for functionObjects that solve PDEs - for finite volume

SourceFiles
    solvePDECommonFiniteVolume.C

Contributors/Copyright:
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef solvePDECommonFiniteVolume_H
#define solvePDECommonFiniteVolume_H

#include "solvePDECommon.H"

namespace Foam
{

// Forward declaration of classes
class FieldValueExpressionDriver;

class solvePDECommonFiniteVolume
:
    public solvePDECommon
{
    // Private Data

        //- The actual field
        autoPtr<volScalarField> theField_;


protected:

    // Protected Data

        //- The driver that calculates the coefficients
        autoPtr<FieldValueExpressionDriver> driver_;


    // Protected Member Functions

        //- Get the actual field
        volScalarField& theField();

        //- Write specific data
        virtual void writeData();

        virtual void writeOldField();

        virtual void writeNewField();


public:

    //- Runtime type information
    TypeName("solvePDECommonFiniteVolume");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        solvePDECommonFiniteVolume
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        solvePDECommonFiniteVolume(const solvePDECommonFiniteVolume&) = delete;


    //- Destructor
    virtual ~solvePDECommonFiniteVolume();


    // Member Functions

        //- Read the forces data
        virtual void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solvePDECommonFiniteVolume&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
