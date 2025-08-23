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

Contributors/Copyright:
    2012-2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "CommonValueExpressionDriver.H"

namespace Foam {


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <typename T>
tmp<T> CommonPluginFunction::evaluate(
    const string &content,
    int &consumed
)
{
    evaluateInternal(content,consumed);

    return tmp<T>(result_.getObjectResult<T>(false));
}

template <typename T>
tmp<Field<T>> CommonPluginFunction::evaluate(
    const string &content,
    int &consumed,
    bool isPoint
)
{
    evaluateInternal(content,consumed);

    if (isPoint!=result_.isPoint()) {
        FatalErrorIn("CommonPluginFunction::evaluate")
            << "Expected 'pointiness' of result " << result_.isPoint()
                << " is not as expected " << isPoint
                << endl
                << exit(FatalError);

    }
    return tmp<Field<T>>(result_.getResult<T>(false));
}

inline CommonValueExpressionDriver &CommonPluginFunction::parentDriver()
{
    return const_cast<CommonValueExpressionDriver &>(parentDriver_);
}

inline const fvMesh &CommonPluginFunction::mesh() const
{
    return parentDriver_.mesh();
}

inline ExpressionResult &CommonPluginFunction::result()
{
    return result_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
