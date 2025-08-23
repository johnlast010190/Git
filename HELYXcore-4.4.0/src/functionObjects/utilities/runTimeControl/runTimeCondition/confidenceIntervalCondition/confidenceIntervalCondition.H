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
    (c) 2016 OpenCFD Ltd.
    (c) 2018 Engys Ltd.
    (c) 2015 OpenFOAM Foundation

Class
    Foam::functionObjects::runTimeControls::confidenceIntervalCondition

Description
    Confidence interval run time condition
    Satisfied when confidence interval drops below a specified range.

SourceFiles
    confidenceIntervalCondition.H
    confidenceIntervalCondition.C

\*---------------------------------------------------------------------------*/

#ifndef confidenceIntervalCondition_H
#define confidenceIntervalCondition_H

#include "runTimeControl/runTimeCondition/runTimeCondition/runTimeCondition.H"
#include "db/functionObjects/writeFile/writeFile.H"
#include "primitives/enums/NamedEnum.H"
#include "primitives/bools/Switch/Switch.H"
#include "primitives/functions/Function2/Function2/Function2.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                Class confidenceIntervalCondition Declaration
\*---------------------------------------------------------------------------*/

class confidenceIntervalCondition
:
    public runTimeCondition,
    public writeFile
{
public:

        //- Interval type options
        enum intvTypes
        {
            itRel,
            itAbs
        };


protected:

    // Protected data

        //- Interval type names
        static const NamedEnum<intvTypes,2> intvTypeNames_;

        //- Function Object Dictionary
        dictionary dict_;

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of entries on which to operate
        word entryName_;

        //- Required confidence - 0.01 -- 0.999
        scalar confidence_;

        //- Termination interval
        const scalar interval_;

        //- interval mode: fraction/absolute
        const intvTypes intervalType_;

        //- Minimum samples
        label minSamples_;

        //- table file name
        fileName tDistFile_;

        //- T distribution table
        autoPtr<Function2<scalar>> tDistribution_;

        //- Distribution Values
        DynamicList<scalar> X_;

        //- Time Values
        DynamicList<scalar> T_;

        //- Gaussian Filter Radius
        label filterRadius_;

        //- Tolerance for the slope of the regressed line
        scalar slopeTolerance_;

        //- Switch to deactivate restart from previous data
        Switch restart_;

        //- Flag to insert or not the smooth function in the transient
        // detection.
        bool smoothValueForTransientDetection_;

        //- Flag to insert or not the smooth function in the convergence assesm.
        bool smoothValueForConvergence_;

        //- Entry in case domain flashing with large time steps occur
        scalar minTransientTime_;

    // Protected Member Functions

        //- function to calculate scalar statistics
        void calc
        (
            bool& satisfied
        );

        List<scalar> lowPassFilter(const List<scalar>& X);

        scalar findGlobalMinimum
        (
            const List<scalar>& X,
            const List<scalar>& t,
            const label& iter
        );

        List<scalar> fractionalFilter(const List<scalar>& X);

        scalar boundLimit(const List<scalar>& X, label i);

        List<scalar> reverseStandardError(const List<scalar>& X, label h);

        scalar meanData(const List<scalar>& list, label j, label h);

        scalar stErr(const List<scalar>& list, label j, label h, scalar mean);

        List<scalar> autoCorrelationFactor(const List<scalar>& X);

        scalar autoCorrelationTime
        (
            const List<scalar>& autoCorrFactor
        );

        scalar simpleLinearRegression
        (
            const List<scalar>& X, const List<scalar>& Y
        );

        void computeEquilibratedUncertainty
        (
            const List<scalar>& X,
            const scalar& aCorrTime,
            scalar& equilibratedAverage,
            scalar& equilibratedUncertainty,
            label& uncorrelatedBatches
        );

        scalar getRegressedDataSlope
        (
            const List<scalar>& T,
            const List<scalar>& X
        );

        bool checkConvergence
        (
            scalar tolerance,
            scalar confidenceInterval,
            scalar slope,
            scalar size
        ) const;

public:

    //- Runtime type information
    TypeName("confidenceInterval");

    //- Constructor
    confidenceIntervalCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~confidenceIntervalCondition();


    // Public Member Functions

        //- Apply the condition
        virtual bool apply(bool postProcess);

        //- Analyze signal in post process mode
        void postProcessMode(bool& satisfied);

        //- Write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimeControl/runTimeCondition/confidenceIntervalCondition/confidenceIntervalConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
