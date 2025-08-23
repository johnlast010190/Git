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
    SolverInfo

Description
    Determines the type of the Solver (whether it is
    compressible/incompressible, turbulent etc)

SourceFiles
    SolverInfo.C

Contributors/Copyright:
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef SolverInfo_H
#define SolverInfo_H

#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class SolverInfo Declaration
\*---------------------------------------------------------------------------*/

class SolverInfo
{
    // Private Data

        word phi_;

        bool compressible_;

        bool turbulent_;

        bool les_;


public:

    // Constructors

        //- Construct from basic information and try to determine properties
        SolverInfo(const objectRegistry &obr,const dictionary& dict);

        //- Disallow default bitwise copy construction
        SolverInfo(const SolverInfo&) = delete;


    // Member Functions

        word phi() const
        {
            return phi_;
        }

        bool compressible() const
        {
            return compressible_;
        }

        bool turbulent() const
        {
            return turbulent_;
        }

        bool les() const
        {
            return les_;
        }

        void printInfo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SolverInfo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
