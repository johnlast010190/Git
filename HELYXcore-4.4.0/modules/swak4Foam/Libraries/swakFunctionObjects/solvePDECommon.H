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
    solvePDECommon

Description
    Common functionality for functionObjects that solve PDEs

SourceFiles
    solvePDECommon.C

Contributors/Copyright:
    2011, 2013-2014, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef solvePDECommon_H
#define solvePDECommon_H

#include "fields/Fields/primitiveFieldsFwd.H"
#include "meshes/primitiveShapes/point/pointFieldFwd.H"
#include "primitives/enums/NamedEnum.H"
#include "fields/volFields/volFields.H"

#include "include/swak.H"

#include "helpers/exprString.H"

#ifdef FOAM_HAS_FVOPTIONS
#ifdef FOAM_FVOPTIONS_IN_FV
#include "cfdTools/general/fvOptions/fvOptionList.H"
#else
#include "fvIOoptionList.H"
#endif
#endif

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyTopoChangeMap;
class polyMeshMap;
class polyMesh;

class solvePDECommon
{
    // Private Data
#ifdef FOAM_HAS_FVOPTIONS

        //- Standin if the solver did not define an option list
        fv::optionList dummyOptionList_;
        bool warnedAboutMissingOptionList_;

#endif

public:

    enum solveAt
    {
        saStartup,
        saTimestep,
        saWrite,
        saNever
    };


protected:

    // Protected Data

        //- Are we working?
        bool active_;

        //- Mesh where the set goes
        const objectRegistry& obr_;

        //- Name of the function object
        const word name_;

        //- Name of the field that is being solved for
        word fieldName_;

        //- Is this a steady calculation
        bool steady_;

        //- Should unsteady runs allow relaxation too
        bool relaxUnsteady_;

        //- The last iteration will also be relaxed
        bool relaxLastIteration_;

        //- If a steady solution does not converge: restore it
        bool restoreNonConvergedSteady_;


    // Protected Member Functions

        //- Should we relax?
        bool doRelax(bool last);

        //- Read expression and dimension
        void readExpressionAndDimension
        (
            const dictionary& dict,
            const word& name,
            exprString& expr,
            dimensionSet& dim
        );

#ifdef FOAM_HAS_FVOPTIONS
        fv::optionList& fvOptions() const;
#endif

        bool needsRhoField(bool warnIfSteady=false) const;


private:

    // Private Member Functions

        //- Do the actual solving
        virtual void solve() = 0;


protected:

    // Protected Data

        static const NamedEnum<solveAt, 4> solveAtNames_;
        solveAt solveAt_;

        //- Should the fields before and after be written?
        bool writeBeforeAfter_;

    // Protected Member Functions

        //- Wrapper around solve()
        void solveWrapper();

        //- write specific data
        virtual void writeData() = 0;

        virtual void writeOldField() = 0;

        virtual void writeNewField() = 0;

        inline bool restoreNonConvergedSteady()
        {
            return steady_ && restoreNonConvergedSteady_;
        }


public:

    //- Runtime type information
    TypeName("solvePDECommon");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        solvePDECommon
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        solvePDECommon(const solvePDECommon&) = delete;


    //- Destructor
    virtual ~solvePDECommon();


    // Member Functions

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Execute
        virtual void execute();

        //- Execute at the final time-loop
        virtual void end();

        //- Write the forces
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
        virtual bool write();
#else
        virtual void write();
#endif

        //- Update for changes of mesh
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update for changes of mesh
        virtual void mapMesh(const polyMeshMap&)
        {}

    //- Update for changes of mesh
#ifdef FOAM_MOVEPOINTS_GETS_POLYMESH
        virtual void movePoints(const polyMesh&)
#else
        virtual void movePoints(const pointField&)
#endif
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solvePDECommon&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
