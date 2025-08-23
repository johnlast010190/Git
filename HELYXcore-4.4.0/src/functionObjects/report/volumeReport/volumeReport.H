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
    Foam::functionObjects::volumeReport

Group
    grpReportFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call volumeReport.write() to calculate&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    volumeReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volumeReport_H
#define functionObjects_volumeReport_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "fieldInstance/fieldInstance.H"
#include "convergenceTermination/convergenceTermination.H"
#include "global/columnatedFileWriter/columnatedFileWriter.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class volumeReport Declaration
\*---------------------------------------------------------------------------*/

class volumeReport: public fvMeshFunctionObject
{
protected:

    // Protected data

        //- input
        Switch logToFile_;

        wordList fields_;
        const dictionary dict_;

        //- Result values
        PackedList<1> reportCells_;
        scalar totalVolume_;
        List< FixedList<scalar, 2>> globalData_;
        List< FixedList<fieldInstance, 2>> localData_;

        //Check for convergence to abort calculation
        PtrList<convergenceTermination> convergenceChecks_;

        columnatedFileWriter fileWriter_;


    // Private Member Functions

        //- Output file header information
        void writeFileHeader() const;

        template<class GeoField>
        void calcNonScalarData
        (
            const scalarField&,
            const volVectorField&,
            const GeoField& f,
            label fI
        );

        void calcData
        (
            const scalarField&,
            const volVectorField&,
            const volScalarField& f,
            label fI
        );

        void updateReportCells();
        void computeTotalVolume();


public:

    //- Runtime type information
    TypeName("volumeReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        volumeReport
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        volumeReport(const volumeReport&) = delete;


    //- Destructor
    virtual ~volumeReport();


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
        void operator=(const volumeReport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
