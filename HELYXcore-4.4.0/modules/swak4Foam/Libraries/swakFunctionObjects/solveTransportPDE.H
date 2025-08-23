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
    solveTransportPDE

Description
    Function object that solves the transport equation

SourceFiles
    solveTransportPDE.C

Contributors/Copyright:
    2011, 2013-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef solveTransportPDE_H
#define solveTransportPDE_H

#include "solvePDECommonFiniteVolume.H"

namespace Foam
{

class solveTransportPDE
:
    public solvePDECommonFiniteVolume
{
protected:

    // Protected Data

        //- The coefficent for the time-derivative
        exprString rhoExpression_;
        dimensionSet rhoDimension_;

        //- The coefficient for the laplacian term
        exprString diffusionExpression_;
        dimensionSet diffusionDimension_;

        //- The name of the phi-field.
        exprString phiExpression_;
        dimensionSet phiDimension_;

        //- The right hand side
        exprString sourceExpression_;
        dimensionSet sourceDimension_;

        //- The right hand side
        exprString sourceImplicitExpression_;
        dimensionSet sourceImplicitDimension_;

        //- Discretize Implicit Source term with SpUp
        bool sourceImplicitUseSuSp_;

        //- Do the actual solving
        void solve();

        //- Should the flux be made relative on a moving mesh
        bool makePhiRelative_;

        //- Velocity for making the flux relative
        exprString velocityExpression_;
        dimensionSet velocityDimension_;
        word velocityName_;


public:

    //- Runtime type information
    TypeName("solveTransportPDE");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        solveTransportPDE
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        solveTransportPDE(const solveTransportPDE&) = delete;


    //- Destructor
    virtual ~solveTransportPDE();


    // Member Functions

        //- Read the forces data
        virtual void read(const dictionary&);
        virtual void setPostProcess(const bool&){};


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solveTransportPDE&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
