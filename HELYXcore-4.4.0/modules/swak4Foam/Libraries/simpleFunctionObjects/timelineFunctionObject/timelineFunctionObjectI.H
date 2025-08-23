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

Contributors/Copyright:
    2008-2011, 2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "timelineFunctionObject.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
/*---------------------------------------------------------------------------*\
                           Class timelineFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<>
inline void timelineFunctionObject::writeData(const word &name,const UList<string> &data) {
    Dbug<< "write Data - string: " << name << " size: " << data.size() << endl;

    OFstream &s=*filePtrs_[name];
    Dbug<< "To: " << s.name() << endl;

    forAll(data,index) {
        string d(data[index]);
        d.replaceAll(" ","_");
        switch(outputFileMode_) {
            case ofmFoam:
            case ofmRaw:
                s << setw
                    (
                        max(
                            label(IOstream::defaultPrecision() + 7),
                            label(d.size()+1)   // at least one space
                        )
                    ) << d.c_str();
                break;
            case ofmCsv:
                s << separator_.c_str() << d.c_str();
                break;
            default:
                FatalErrorIn("timelineFunctionObject::write(const word &name,const UList<T> &data)")
                    << "Don't know how to write files of type "
                        << outputFileModeNames_[outputFileMode_] << endl
                        << exit(FatalError);
        }
    }
}

template<>
inline void timelineFunctionObject::writeData(const word &name,const UList<scalar> &data) {
    Dbug<< "write Data - scalar: " << name << " size: " << data.size() << endl;

    OFstream &s=*filePtrs_[name];
    Dbug<< "To: " << s.name() << endl;

    forAll(data,index) {
        switch(outputFileMode_) {
            case ofmFoam:
            case ofmRaw:
                s << setw(IOstream::defaultPrecision() + 7) << data[index];
                break;
            case ofmCsv:
                s << separator_.c_str() << data[index];
                break;
            default:
                FatalErrorIn("timelineFunctionObject::write(const word &name,const UList<T> &data)")
                    << "Don't know how to write files of type "
                        << outputFileModeNames_[outputFileMode_] << endl
                        << exit(FatalError);
        }
    }
}

template<class T>
void timelineFunctionObject::writeData(const word &name,const UList<T> &data) {
    Dbug<< "write Data: " << name << " type: " << pTraits<T>::typeName
        << " size: " << data.size() << endl;

    OFstream &s=*filePtrs_[name];
    Dbug<< "To: " << s.name() << endl;

    forAll(data,index) {
        if (outputFileMode_==ofmFoam) {
            s << setw(IOstream::defaultPrecision() + 7) << data[index];
        } else {
            for (label i=0;i<pTraits<T>::nComponents;i++) {
                switch(outputFileMode_) {
                    case ofmRaw:
                        s << setw(IOstream::defaultPrecision() + 7) << data[index][i];
                        break;
                    case ofmCsv:
                        s << separator_.c_str() << data[index][i];
                        break;
                    default:
                        FatalErrorIn("timelineFunctionObject::write(const word &name,const UList<T> &data)")
                            << "Don't know how to write files of type "
                                << outputFileModeNames_[outputFileMode_] << endl
                                << exit(FatalError);
                }
            }
        }
    }
}

template<class T>
void timelineFunctionObject::writeTime(const word &name,const T &time) {
    Dbug<< "writeTime: " << name << " t=" << time << endl;
    OFstream &s=*filePtrs_[name];
    Dbug<< "To: " << s.name() << endl;

    switch(outputFileMode_) {
        case ofmFoam:
        case ofmRaw:
            s << setw(IOstream::defaultPrecision() + 7) << time;
            break;
        case ofmCsv:
            s <<  time;
            break;
        default:
            FatalErrorIn("timelineFunctionObject::writeTime(const word &name,const UList<T> &data)")
                    << "Don't know how to write files of type "
                        << outputFileModeNames_[outputFileMode_] << endl
                        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
