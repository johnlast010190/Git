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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::functionObjects::fieldAverage

Group
    grpFieldFunctionObjects

Description
    Calculates average quantities for a user-specified selection of volumetric
    and surface fields.
    With the %subRegion option, also supports fields on a surfMesh.

    Fields are entered as a list of sub-dictionaries, which indicate the type of
    averages to perform, and can be updated during the calculation. The current
    options include:
    - \c mean: arithmetic mean
        \f[
            \overline{x} = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N x_i
        \f]
    - \c prime2Mean: prime-squared mean
        \f[
            \overline{x'}^2 = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N
            (x_i - \overline{x})^2
        \f]
    - \c base: average over 'time', or 'iteration' (\f$N\f$ in the above)
    - \c window: optional averaging window, specified in 'base' units

    Average field names are constructed by concatenating the base field with
    the averaging type, e.g. when averaging field 'U', the resultant fields
    are:
    - arithmetic mean field, \c UMean
    - prime-squared field, \c UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a per-field basis to the
    \c "<functionObject name>Properties" dictionary, located in
    \<time\>/uniform.

    When restarting form a previous calculation, the averaging is continuous or
    may be restarted using the \c restartOnRestart option.

    The averaging process may be restarted after each calculation output time
    using the \c restartOnOutput option or restarted periodically using the \c
    periodicRestart option and setting \c restartPeriod to the required
    averaging period.

    When 'timeIntegral' is enabled for a specific field, the mean is replaced by
    the time average of that field (see example in the usage).

Usage
    \verbatim
    fieldAverage1
    {
        type                fieldAverage;
        libs                ("libfieldFunctionObjects.so");

        writeControl        writeTime;

        restartOnRestart    false;
        restartOnOutput     false;
        periodicRestart     false;
        restartPeriod       0.002;

        fields
        (
            U
            {
                mean            on;
                prime2Mean      on;
                base            time;
                window          10.0;
                windowName      w1;
            }
            p
            {
                mean            on;
                prime2Mean      on;
                base            time;
            }
            alpha.water
            {
                mean            on;
                prime2Mean      off;
                base            time;
                timeIntegral    true; // compute time integral instead of mean
            }
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property          | Description               | Required | Default
        type              | type name: fieldAverage              | yes |
        restartOnRestart  | Restart the averaging on restart     | no  | no
        restartOnOutput   | Restart the averaging on output      | no  | no
        periodicRestart   | Periodically restart the averaging   | no  | no
        restartPeriod     | Periodic restart period              | conditional |
        restartTime       | One-shot reset of the averaging | no | great
        fields            | list of fields and averaging options | yes |
        subRegion         | name of a sub-region such as a surface name | no |
        Welford           | Welford formula for prime2Mean       | no | no
        symmetrize        | symmetrize the Prime2Mean of vectors | no | yes
    \endtable

Notes
    - To employ the \c prime2Mean option, the \c mean option must be selected.
    - Suggested use of Welford formula for prime2Mean for reduced round-off
      error especially when running in single-precision (SP).
    - Symmetrize option is only valid when Welford is enabled.

References
    [1] B.P Welford:
        "Note on a method for calculating corrected sums of squares and
        products", Technometrics. 4 (3): 419–420. doi:10.2307/1266577. 1962.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    fieldAverage.C
    fieldAverageTemplates.C
    fieldAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldAverage_H
#define functionObjects_fieldAverage_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fieldAverageItem/fieldAverageItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldAverage
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- Restart the averaging process on restart
        Switch restartOnRestart_;

        //- Restart the averaging process on output
        Switch restartOnOutput_;

        //- Periodically restart the averaging process
        Switch periodicRestart_;

        //- Restart period
        scalar restartPeriod_;

        //- Specific restart time
        scalar restartTime_;

        //- Welford formula
        Switch Welford_;

        //- Symmetrize the prime2Mean if Welford formula is used.
        //  Prime2Mean of vector in Welford formula results in tensor.
        //  Use this flag to symmetrize the tensor. Off diagonal elements will
        //  be almost zero.
        Switch symmetrize_;

        //- List of field average items, describing what averages to be
        //  calculated and output.
        List<fieldAverageItem> faItems_;

        //- Iteration steps counter
        List<label> totalIter_;

        //- Total time counter
        List<scalar> totalTime_;

        //- Index for periodic restart
        label periodIndex_;


    // Protected Member Functions

        // Initialisation routines

            //- Reset lists (clear existing values) and initialise averaging.
            //  Check requested field averages are valid, populate field lists.
            void initialise();

            //- Read the mean average field
            template<class Type>
            void readMeanFieldType(const label fieldi);

            //- Read the mean average field
            template<class Type>
            void readMeanField(const label fieldi);

            //- Initialise the mean average field
            template<class Type>
            void initialiseMeanFieldType(const label fieldi);

            //- Initialise the mean average field
            template<class Type>
            void initialiseMeanField(const label fieldi);

            //- Read the prime-squared average field
            template<class Type1, class Type2>
            void readPrime2MeanFieldType(const label fieldi);

            //- Read the prime-squared average field
            template<class Type1, class Type2>
            void readPrime2MeanField(const label fieldi);

            //- Initialise the prime-squared average field
            template<class Type1, class Type2>
            void initialisePrime2MeanFieldType(const label fieldi);

            //- Initialise the prime-squared average field
            template<class Type1, class Type2>
            void initialisePrime2MeanField(const label fieldi);

            // Welford templated functions.
            // Calculate prime-squared average fields.

                //- Initialise the Welford prime-squared average field
                template<class Type1, class Type2>
                void initialisePrime2MeanFieldWelType(const label fieldi);

                //- Initialise the Welford prime-squared average field
                template<class Type1, class Type2>
                void initialisePrime2MeanFieldWel(const label fieldi);

                //- Initialise the vector-symmTensor field
                template<class Type1, class Type2>
                void initialisePrime2MeanFieldWelSymmType(const label fieldi);

                //- Initialise the vector-symmTensor field
                template<class Type1, class Type2>
                void initialisePrime2MeanFieldWelSymm(const label fieldi);

            //- Add duplicate type field to database
            template<class Type>
            void copyFieldType(const word&, const word&, bool&, bool&);

            //- Add min average field to database
            template<class Type>
            void addMinField(const label fieldi);

            //- Add max average field to database
            template<class Type>
            void addMaxField(const label fieldi);


        // Calculation functions

            //- Main calculation routine
            void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFieldType(const label fieldi) const;

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields() const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFieldType(const label fieldi) const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields() const;

            // Welford templated functions

                //- Calculate prime-squared average fields
                template<class Type1, class Type2>
                void calculateMeanPrime2MeanFieldWelType(const label) const;

                //- Instantiate for vector-symmTensor
                template<class Type1, class Type2>
                void calculateMeanPrime2MeanFieldWelSymmType(const label) const;

                //- Calculate prime-squared average fields
                template<class Type1, class Type2>
                void calculateMeanPrime2MeanFieldsWel() const;

                //- Instantiate for vector-symmTensor
                template<class Type1, class Type2>
                void calculateMeanPrime2MeanFieldsWelSymm() const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2MeanType(const label fieldi) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean() const;

            //- Calculate min average fields
            template<class Type>
            void calculateExtFieldType
            (
                const word& fieldName,
                const word& outputFieldName,
                scalar
            ) const;

            //- Calculate min average fields
            template<class Type>
            void calculateMinFields() const;

            //- Calculate Max average fields
            template<class Type>
            void calculateMaxFields() const;


        // I-O

            //- Write averages
            void writeAverages();

            //- Write fields
            template<class Type>
            void writeFieldType(const word& fieldName) const;

            //- Write fields
            template<class Type>
            void writeFields() const;

            //- Read
            void read(const dictionary& dict, const bool construct);


public:

    //- Runtime type information
    TypeName("fieldAverage");


    // Constructors

        //- Construct from name, Time and dictionary
        fieldAverage(const word& name, const Time& runTime,const dictionary&);

        //- Disallow default bitwise copy construction
        fieldAverage(const fieldAverage&) = delete;


    //- Destructor
    virtual ~fieldAverage();


    // Member Functions

        //- Restart averaging for restartOnOutput
        void restart();

        //- Read the field average data
        virtual bool read(const dictionary&) override;

        //- Calculate the field averages
        virtual bool execute() override;

        //- Write the field averages
        virtual bool write() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldAverage&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldAverage/fieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
