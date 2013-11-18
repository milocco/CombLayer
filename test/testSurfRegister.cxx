/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   test/testSurfRegister.cxx
*
 * Copyright (c) 2004-2013 by Stuart Ansell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Triple.h"
#include "Surface.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "surfIndex.h"
#include "surfRegister.h"

#include "testFunc.h"
#include "testSurfRegister.h"

using namespace ModelSupport;

testSurfRegister::testSurfRegister()
  /*!
    Constructor
  */
{}

testSurfRegister::~testSurfRegister() 
  /*!
    Destructor
  */
{}


int 
testSurfRegister::applyTest(const int extra)
  /*!
    Applies all the tests and returns 
    the error number
    \param extra :: Index of test
    \returns -ve on error 0 on success.
  */
{
  ELog::RegMethod RegA("testSurfRegister","applyTest");
  TestFunc::regSector("testSurfRegister");

  typedef int (testSurfRegister::*testPtr)();
  testPtr TPtr[]=
    {
      &testSurfRegister::testIdentical,
      &testSurfRegister::testPlaneReflection,
      &testSurfRegister::testUnique
    };

  const std::string TestName[]=
    {
      "Identical",
      "PlaneReflection",
      "Unique"
    };

  const int TSize(sizeof(TPtr)/sizeof(testPtr));
  if (!extra)
    {
      std::ios::fmtflags flagIO=std::cout.setf(std::ios::left);
      for(int i=0;i<TSize;i++)
        {
	  std::cout<<std::setw(30)<<TestName[i]<<"("<<i+1<<")"<<std::endl;
	}
      std::cout.flags(flagIO);
      return 0;
    }
  for(int i=0;i<TSize;i++)
    {
      if (extra<0 || extra==i+1)
        {
	  TestFunc::regTest(TestName[i]);
	  const int retValue= (this->*TPtr[i])();
	  if (retValue || extra>0)
	    return retValue;
	}
    }
  return 0;
}

int
testSurfRegister::testIdentical()
  /*!
    Test to see if a number of equal surfaces are all assigned 
    to the correct number.
    \returns 0 on succes and -ve on failure
  */
{
  ELog::RegMethod RegA("testSurfRegister","testIdentical");

  // Start from a clean sheet
  surfRegister SMap;
  ModelSupport::surfIndex& SurI=ModelSupport::surfIndex::Instance();
  SurI.reset();

  typedef boost::tuple<int,std::string,int> TTYPE;

  std::vector<TTYPE> Tests;
  Tests.push_back(TTYPE(3,"px 34",3));
  Tests.push_back(TTYPE(4,"px 34",3));
  Tests.push_back(TTYPE(5,"px 34",3));
  Tests.push_back(TTYPE(6,"px 34",3));

    
  std::vector<TTYPE>::const_iterator tc;
  int flag(0);
  for(tc=Tests.begin();tc!=Tests.end() && !flag;tc++)
    {
      const int SN(tc->get<0>());
      Geometry::Plane* PX=SurI.createUniqSurf<Geometry::Plane>(SN);

      if (!PX || PX->setSurface(tc->get<1>()))
	{
	  ELog::EM<<"Failed : surf == "<<tc->get<0>()<<ELog::endCrit;
	  if (PX) 
	    ELog::EM<<"PX      "<<*PX<<ELog::endCrit;
	  else
	    ELog::EM<<"PX NULL "<<ELog::endCrit;
	  return -1;
	}

      if (SMap.registerSurf(SN,PX)!=tc->get<2>())
	{
	  ELog::EM<<"Surface number == "<<SMap.realSurf(SN)<<ELog::endCrit;
	  flag=-2;
	}
    }
  return flag;
}

int
testSurfRegister::testPlaneReflection()
  /*!
    Test to see if a number of equal surfaces are all assigned 
    to the correct number.
    \returns 0 on succes and -ve on failure
  */
{
  ELog::RegMethod RegA("testSurfRegister","testPlaneReflection");

  // Start from a clean sheet
  surfRegister SMap;
  ModelSupport::surfIndex& SurI=ModelSupport::surfIndex::Instance();
  SurI.reset();

  typedef boost::tuple<int,std::string,int> TTYPE;

  std::vector<TTYPE> Tests;
  Tests.push_back(TTYPE(3,"px 1",3));
  Tests.push_back(TTYPE(4,"p -1 0 0 -1",3));

    
  std::vector<TTYPE>::const_iterator tc;
  int flag(0);
  for(tc=Tests.begin();tc!=Tests.end() && !flag;tc++)
    {
      const int SN(tc->get<0>());
      Geometry::Plane* PX=SurI.createUniqSurf<Geometry::Plane>(SN);

      if (!PX || PX->setSurface(tc->get<1>()))
	{
	  ELog::EM<<"Failed : surf == "<<tc->get<0>()<<ELog::endCrit;
	  if (PX) 
	    ELog::EM<<"PX      "<<*PX<<ELog::endCrit;
	  else
	    ELog::EM<<"PX NULL "<<ELog::endCrit;
	  return -1;
	}

      if (SMap.registerSurf(SN,PX)!=tc->get<2>())
	{
	  ELog::EM<<"Surface number == "<<SMap.realSurf(SN)<<ELog::endCrit;
	  Geometry::Plane* PY=SMap.realPtr<Geometry::Plane>(tc->get<2>());
	  if (PY)
	    ELog::EM<<"Plane[Expected] == "<<*PY<<ELog::endDebug;
	  ELog::EM<<"Plane == "<<*PX<<ELog::endDebug;
	  ELog::EM<<"Surf == "<<PX->side(Geometry::Vec3D(1,0,0))
		  <<ELog::endDebug;
	  ELog::EM<<"Surf == "<<PY->side(Geometry::Vec3D(1,0,0))
		  <<ELog::endDebug;
	  PX->mirrorSelf();
	  ELog::EM<<"Rev Surf == "<<PX->side(Geometry::Vec3D(1,0,0))
		  <<ELog::endDebug;
	  flag=-2;
	}
    }
  return flag;
}

int
testSurfRegister::testUnique()
  /*!
    Test to see if a number of equal surfaces are all assigned 
    to the correct number.
    \returns 0 on succes and -ve on failure
  */
{
  ELog::RegMethod RegA("testSurfRegister","testUnique");

  // Start from a clean sheet
  surfRegister SMap;
  ModelSupport::surfIndex& SurI=ModelSupport::surfIndex::Instance();
  SurI.reset();

  typedef boost::tuple<int,std::string,int> TTYPE;
  typedef boost::tuple<int,double> RTYPE;   // Result checkx

  std::vector<TTYPE> Tests;
  Tests.push_back(TTYPE(3,"px 34",3));
  Tests.push_back(TTYPE(4,"px 35",4));
  Tests.push_back(TTYPE(5,"px 36",5));
  Tests.push_back(TTYPE(6,"px 34",3));
  Tests.push_back(TTYPE(6,"px 34",3));
  // Result check
  std::vector<RTYPE> ResTest;
  ResTest.push_back(RTYPE(3,34));
  ResTest.push_back(RTYPE(4,35));
  ResTest.push_back(RTYPE(5,36));
  ResTest.push_back(RTYPE(6,34));

    
  std::vector<TTYPE>::const_iterator tc;
  int flag(0);
  for(tc=Tests.begin();tc!=Tests.end() && !flag;tc++)
    {
      const int SN(tc->get<0>());
      Geometry::Plane* PX=SurI.createUniqSurf<Geometry::Plane>(SN);
      if (!PX || PX->setSurface(tc->get<1>()))
	{
	  ELog::EM<<"Failed : surf == "<<tc->get<0>()<<ELog::endCrit;
	  if (PX) 
	    ELog::EM<<"PX      "<<*PX<<ELog::endCrit;
	  else
	    ELog::EM<<"PX NULL "<<ELog::endCrit;
	  return -1;
	}

      if (SMap.registerSurf(SN,PX)!=tc->get<2>())
	{
	  ELog::EM<<"Surface number == "<<SMap.realSurf(SN)<<ELog::endCrit;
	  return -2;
	}
    }
  std::vector<RTYPE>::const_iterator rc;
  for(rc=ResTest.begin();rc!=ResTest.end();rc++)
    {
      const Geometry::Plane* PN=
	dynamic_cast<const Geometry::Plane*>(SMap.realSurfPtr(rc->get<0>()));
      if (!PN || fabs(PN->getDistance()-rc->get<1>())>1e-7)
	{
	  ELog::EM<<"Failed Surface number == "<<rc->get<0>()<<ELog::endCrit;
	  ELog::EM<<"Expect == "<<rc->get<1>()<<ELog::endCrit;
	  if (PN)
	    ELog::EM<<"PN == "<<*PN<<ELog::endCrit;
	  return -3;
	}
    }
  return 0;
}

