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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::LagrangianCloudSourcePluginFunction

Description

SourceFiles
    LagrangianCloudSourcePluginFunction.C

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef LagrangianCloudSourcePluginFunction_H
#define LagrangianCloudSourcePluginFunction_H

#include "plugins/FieldValuePluginFunction.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class LagrangianCloudSourcePluginFunction Declaration
\*---------------------------------------------------------------------------*/

class LagrangianCloudSourcePluginFunction
:
    public FieldValuePluginFunction
{
    // Private Data

        word cloudName_;
        HashSet<word> supportedTypes_;


protected:

    // Protected Member Functions

        template<class LookupCloudType>
        bool hasCloudAs();

        Ostream& listAvailableClouds(Ostream& o);

        template<class OtherCloudType, class LookupCloudType>
        const OtherCloudType& getCloudAs();

        const word& cloudName() const
        {
            return cloudName_;
        }

        const HashSet<word>& supportedTypes() const
        {
            return supportedTypes_;
        }


public:

    // Constructors

        //- Construct from components
        LagrangianCloudSourcePluginFunction
        (
            const FieldValueExpressionDriver& parentDriver,
            const word& name,
            const word& returnValueType,
            const string& additionalArgs = ""
        );

        //- Disallow default bitwise copy construction
        LagrangianCloudSourcePluginFunction
        (
            const LagrangianCloudSourcePluginFunction&
        ) = delete;


    //- Destructor
    virtual ~LagrangianCloudSourcePluginFunction()
    {}


    // Member Functions

        void setArgument(label index, const word&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianCloudSourcePluginFunction&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LagrangianCloudSourcePluginFunctionI.H"
#endif

// return a autoPtr to the needed field. This is required because
// "similar" clouds don't really inherit (because of the template paramter)
#define tryCall(rType,neededType,usedType,theCall)    \
    if (hasCloudAs<neededType>()) {                             \
        return autoPtr<rType> (                                \
            new rType(                                         \
                getCloudAs<neededType,usedType>().theCall      \
            )                                                  \
        );                                                     \
    }

// Fail if no cloud data was found
#define noCloudFound(pVar) \
    if (!pVar.valid()) {                                                 \
        listAvailableClouds(Info);                                      \
        FatalErrorIn("Plugin function "+this->name())                   \
            << "No matching cloud " << this->cloudName() << " found"    \
                << endl                                                 \
                << "Supported cloud types (typeid): "                   \
                << supportedTypes() << endl                             \
                << exit(FatalError);                                    \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
