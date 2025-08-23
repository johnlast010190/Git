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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::noiseModel

Description
    Base class for noise models.

    Data is read from a dictionary, e.g.

    \verbatim
    rhoRef          1;
    N               4096;
    fl              25;
    fu              25;
    startTime       0;

    outputPrefix    "test1";

    // Optional write options dictionary
    writeOptions
    {
        writePrmsf  no;
        writeSPL    yes;
        writePSD    yes;
        writePSDf   no;
        writeOctaves yes;
    }
    \endverbatim

    where
    \table
        Property     | Description                   | Required  | Default value
        rhoRef       | Reference density             | no        | 1
        N            | Number of samples in sampling window | no | 65536 (2^16)
        fl           | Lower frequency bounds        | no        | 25
        fu           | Upper frequency bounds        | no        | 10000
        startTime    | Start time                    | no        | 0
        pType        | type of noise pressure        | no        | static [dynamic, p/rho]
        p0           | reference sound pressure      | no        | 2e-5
        outputPrefix | Prefix applied to output files| no        | ''
        graphFormat  | Graph format                  | no        | raw
        writePrmsf   | Write Prmsf data              | no        | yes
        writeSPL     | Write SPL data                | no        | yes
        writePSD     | Write PSD data                | no        | yes
        writePSDf    | Write PSDf data               | no        | yes
        writeOctaves | Write octaves data            | no        | yes
    \endtable

Note
    The number of samples in the sampling window must be a power of 2


SourceFiles
    noiseModel.C

\*---------------------------------------------------------------------------*/

#ifndef noiseModel_H
#define noiseModel_H

#include "db/dictionary/dictionary.H"
#include "primitives/Scalar/lists/scalarList.H"
#include "db/Time/instant/instantList.H"
#include "windowModels/windowModel/windowModel.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noiseModel Declaration
\*---------------------------------------------------------------------------*/

class noiseModel
{
protected:

    // Protected Data

        //- Copy of dictionary used for construction
        const dictionary dict_;

        //- Reference density (to convert from kinematic to static pressure)
        scalar rhoRef_;

        //- Number of samples in sampling window, default = 2^16
        label nSamples_;

        //- Lower frequency limit, default = 25Hz
        scalar fLower_;

        //- Upper frequency limit, default = 10kHz
        scalar fUpper_;

        //- Flag to indicate that custom frequency bounds are being used
        bool customBounds_;

        //- Start time, default = 0s
        scalar startTime_;

        //- Window model
        autoPtr<windowModel> windowModelPtr_;

        //- Graph format
        word graphFormat_;


        // Data validation

            //- Min pressure value
            scalar minPressure_;

            //- Min pressure value
            scalar maxPressure_;

			//pressure type
			word pType_;

			// reference sound pressure
			scalar p0_;

        // Write options

            //- Output  file prefix, default = ''
            fileName outputPrefix_;

            //- Write Prmsf; default = yes
            bool writePrmsf_;

            //- Write SPL; default = yes
            bool writeSPL_;

            //- Write PSD; default = yes
            bool writePSD_;

            //- Write PSDf; default = yes
            bool writePSDf_;

            //- Write writeOctaves; default = yes
            bool writeOctaves_;


    // Protected Member Functions

        //- Helper function to read write options and provide info feedback
        void readWriteOption
        (
            const dictionary& dict,
            const word& lookup,
            bool& option
        ) const;

        //- Check and return uniform time step
        scalar checkUniformTimeStep
        (
            const scalarList& times
        ) const;

		scalar checkUniformTimeStep
        (
            const scalarList& times,
            bool &isUniform
        ) const;

        //- Return true if all pressure data is within min/max bounds
        bool validateBounds(const scalarList& p) const;

        //- Find and return start time index
        label findStartTimeIndex
        (
            const instantList& allTimes,
            const scalar startTime
        ) const;

        //- Return the base output directory
        fileName baseFileDir(const label dataseti) const;


public:

    //- Runtime type information
    TypeName("noiseModel");

    //- Run time selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        noiseModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Selector
    static autoPtr<noiseModel> New(const dictionary& dict);

    // Constructors

        //- Construct from components
        noiseModel(const dictionary& dict, const bool readFields = true);

        //- Disallow default bitwise copy construction
        noiseModel(const noiseModel&) = delete;


    //- Destructor
    virtual ~noiseModel();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Abstract call to calculate
        virtual void calculate() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noiseModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
