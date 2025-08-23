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
    EvolveCloudFunctionObject

Description
    Basis for the other function objects here. Evaluates the most common options

SourceFiles
    EvolveCloudFunctionObject.C

Contributors/Copyright:
    2012-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#ifndef EvolveCloudFunctionObject_H
#define EvolveCloudFunctionObject_H

#include "db/functionObjects/functionObject/functionObject.H"
#include "db/dictionary/dictionary.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/objectRegistry/objectRegistry.H"
#include "fields/volFields/volFields.H"
#include "include/swak.H"
#include "helpers/switchableNotImplemented.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class EvolveCloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EvolveCloudFunctionObject
:
    public functionObject
{
    // Private Data

        autoPtr<CloudType> cloud_;

        //- Input dictionary
        dictionary dict_;

        //- name of the region
        word regionName_;

        //- the registry
        const objectRegistry& obr_;

        //- name of the cloud
        word cloudName_;

        //- direction of gravity
        dimensionedVector g_;

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        //- Make sure that execute is not called twice at one time-step
        label lastTimeStepExecute_;
#endif

protected:

    // Protected Member Functions

        //- Get a reference to the field
        template <class FieldType >
        const FieldType& getField(const word& fieldName);

        const word& cloudName() const
        {
            return cloudName_;
        }

        autoPtr<CloudType>& cloud()
        {
            return cloud_;
        }

        const dimensionedVector& g() const
        {
            return g_;
        }

        const objectRegistry& obr() const
        {
            return obr_;
        }

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        virtual bool start() = 0;
#endif

public:

    // Constructors

        //- Construct from components
        EvolveCloudFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        EvolveCloudFunctionObject(const EvolveCloudFunctionObject&) = delete;


    // Member Functions

        //- Execute is called at each ++ or += of the time-loop
        virtual bool execute(bool forceWrite);

#ifdef FOAM_FUNCTIONOBJECT_EXECUTE_HAS_NO_FORCE
        bool execute()
        {
            return execute(false);
        }
#endif

#ifdef FOAM_FUNCTIONOBJECT_HAS_SEPARATE_WRITE_METHOD_AND_NO_START
        virtual bool execute()
        {
            return execute(false);
        }

        virtual bool write()
        {
            return execute(true);
        }
#endif

        //- Read and set the function object if it's data has changed
        virtual bool read(const dictionary& dict);

#ifdef FOAM_FUNCTIONOBJECT_NEEDS_UPDATEMESH

        //- Make functionObject happy but fail because
        //  this is not an implementation
        virtual void topoChange(const polyTopoChangeMap&)
        {
            switchableNotImplemented
            (
                "EvolveCloudFunctionObject::topoChange",
                dict_
            );
        }

        //- Make functionObject happy but fail because
        //  this is not an implementation
        virtual void movePoints(const polyMesh&)
        {
            switchableNotImplemented
            (
                "EvolveCloudFunctionObject::movePoints",
                dict_
            );
        }
#endif


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const EvolveCloudFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "EvolveCloudFunctionObject.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
