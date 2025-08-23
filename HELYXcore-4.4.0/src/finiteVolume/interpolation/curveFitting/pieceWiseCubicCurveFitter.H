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
    (c) 2010-2012 Engys Ltd.
 \*---------------------------------------------------------------------------*/
#ifndef _PIECEWISECUBICCURVEFITTER_
#define _PIECEWISECUBICCURVEFITTER_

#include "fields/Fields/Field/Field.H"
#include "fields/Fields/vector2DField/vector2DFieldFwd.H"

namespace Foam
{
  class pieceWiseCubicCurveFitter{

  private :
    vector2DField originalData_;
    /*
    void setInterval( const scalar & input );
    void setCoefficients( );
    */

  public :
    pieceWiseCubicCurveFitter(const vector2DField & data);
    pieceWiseCubicCurveFitter(const pieceWiseCubicCurveFitter & orig){
        this->originalData_ = orig.originalData_;
    };
    virtual ~pieceWiseCubicCurveFitter(){ };

    scalar interpolate( const scalar & input) const;
  };


}

#endif //_PIECEWISECUBICCURVEFITTER_
