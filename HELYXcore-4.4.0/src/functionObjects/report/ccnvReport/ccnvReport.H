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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2024 Engys Ltd.

Class
    Foam::functionObjects::ccnvReport

Group
    grpReportFunctionObjects

Description

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    ccnvReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ccnvReport_H
#define functionObjects_ccnvReport_H

#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include <iostream>
#include <fstream>
#include "containers/Lists/PackedList/PackedList.H"
#include "fieldInstance/fieldInstance.H"
#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class ccnvReport Declaration
\*---------------------------------------------------------------------------*/

class ccnvReport
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- CCNV precision (4 single and 8 double)
        long ccnv_pr_;

        //- Velocity field name
        word UName_;

        //- Phi field name
        word PhiName_;

        //- Pressure field name
        word PName_;

        //- Step number
        int stepNumber_;

        //- Input
        Switch logToFile_;

        //- Log file name
        fileName logFileName_;

        //- Log file ptr
        std::ofstream logFilePtr_;

        //- Result values
        double maxDivergence_;


public:

    //- Runtime type information
    TypeName("ccnvReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        ccnvReport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ccnvReport(const ccnvReport&) = delete;


    //- Destructor
    virtual ~ccnvReport();


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ccnvReport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
