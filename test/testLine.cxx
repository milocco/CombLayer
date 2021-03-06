/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   test/testLine.cxx
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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <list>
#include <vector>
#include <map>
#include <complex>
#include <string>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "mathSupport.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Transform.h"
#include "Surface.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "General.h"
#include "Line.h"
#include "LineIntersectVisit.h"
#include "SurInter.h"

#include "testFunc.h"
#include "testLine.h"

using namespace Geometry;

const double LTolerance(1e-8);

testLine::testLine() 
  /// Constructor
{}

testLine::~testLine() 
  /// Destructor
{}

int 
testLine::applyTest(const int extra)
  /*!
    Applies all the tests and returns 
    the error number
    \param extra :: index of test
    \retval -1 Failed
    \retval 0 All succeeded
  */
{
  ELog::RegMethod RegA("testLine","applyTest");
  TestFunc::regSector("testLine");

  typedef int (testLine::*testPtr)();
  testPtr TPtr[]=
    {
      &testLine::testIntersection,
      &testLine::testInterDistance
    };
  const std::string TestName[]=
    {
      "Intersection",
      "InterDistance"
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
testLine::testConeIntersect()
  /*!
    Test the intersection of cones
    \return -ve on error
  */
{
  ELog::RegMethod RegA("testLine","testCone");

  // Cone : Start Point : Normal : NResults : distance A : distance B 
  typedef boost::tuple<std::string,Geometry::Vec3D,Geometry::Vec3D,
		       size_t,double,double> TTYPE;
  std::vector<TTYPE> Tests;
  
  Tests.push_back(TTYPE("ky 1 1 1",
			Geometry::Vec3D(-3,0,0),Geometry::Vec3D(1,0,0),
			0,0.0,0.0));
  Tests.push_back(TTYPE("ky 1 1 -1",
			Geometry::Vec3D(-3,0,0),Geometry::Vec3D(1,0,0),
			2,4.0,2.0));
  Tests.push_back(TTYPE("ky 12.3 1 -1",
			Geometry::Vec3D(-7.0,8.90,0),Geometry::Vec3D(1,0,0),
			2,10.4,3.6));
  
  std::vector<TTYPE>::const_iterator tc;
  for(tc=Tests.begin();tc!=Tests.end();tc++)
    {
      Cone A;
      Line LX;
      LX.setLine(tc->get<1>(),tc->get<2>());
      
      const int retVal=A.setSurface(tc->get<0>());
      if (retVal)
        {
	  ELog::EM<<"Failed to build "<<tc->get<0>()
		  <<" Ecode == "<<retVal<<ELog::endErr;
	  return -1;
	}
      std::vector<Geometry::Vec3D> OutPt;
      const size_t NR=LX.intersect(OutPt,A);

      if (NR!=tc->get<3>())
	{
	  ELog::EM<<"Failure for test "<<tc-Tests.begin()<<ELog::endCrit;
	  ELog::EM<<"Solution Count"<<NR<<" ["<<
	    tc->get<3>()<<"] "<<ELog::endCrit;
	  return -1;
	}
      const double DA=(NR>0) ? OutPt[0].Distance(tc->get<1>()) : 0.0;
      const double DB=(NR>1) ? OutPt[1].Distance(tc->get<1>()) : 0.0;
      if (NR>0 && fabs(DA-tc->get<4>())> 1e-5) 
	{
	  ELog::EM<<"Failure for test "<<tc-Tests.begin()<<ELog::endCrit;
	  ELog::EM<<"Point A "<<OutPt[0]<<" :: "<<
	    tc->get<1>()+tc->get<2>()*tc->get<4>()<<ELog::endDebug;
	  ELog::EM<<"DA "<<DA<<ELog::endCrit;
	  return -1;
	}
      if (NR>1 && fabs(DB-tc->get<5>())> 1e-5) 
	{
	  ELog::EM<<"Failure for test "<<tc-Tests.begin()<<ELog::endCrit;
	  ELog::EM<<"Point B "<<OutPt[1]<<" :: "<<
	    tc->get<1>()+tc->get<2>()*tc->get<5>()<<ELog::endCrit;
	  ELog::EM<<"DB "<<DB<<ELog::endCrit;
	  return -1;
	}
    }
  return 0;
}

int
testLine::testIntersection()
  /*!
    Tests the Creation of Line and stuff
    \return 0 sucess / -ve on failure
  */
{
  ELog::RegMethod RegA("testLine","testIntersection");
  if (testConeIntersect())
    return -1;
  return 0;
}

int
testLine::testInterDistance()
  /*!
    Test the interPoint/interDistance template function
    \return 0 sucess / -ve on failure
  */
{
  ELog::RegMethod RegItem("testLine","testInterDistance");


  std::vector<boost::shared_ptr<Geometry::Surface> > SurList;
  SurList.push_back(boost::shared_ptr<Geometry::Surface>(new Geometry::Plane(1,0)));
  SurList.back()->setSurface("px 80");
  SurList.push_back(boost::shared_ptr<Geometry::Surface>(new Geometry::Cylinder(2,0)));
  SurList.back()->setSurface("c/z 3 5 50");

  // surfN : Origin : Axis : results 
  typedef boost::tuple<const Geometry::Surface*,Geometry::Vec3D,
		       Geometry::Vec3D,double> TTYPE;
  std::vector<TTYPE> Tests;
  
  Tests.push_back(TTYPE(SurList[0].get(),Geometry::Vec3D(0,0,0),
			Geometry::Vec3D(1,0,0),80.0));

  Tests.push_back(TTYPE(SurList[1].get(),Geometry::Vec3D(0,0,0),
			Geometry::Vec3D(1,0,0),sqrt(50*50-25)+3.0));

  // DO Tests:
  std::vector<TTYPE>::const_iterator tc;
  for(tc=Tests.begin();tc!=Tests.end();tc++)
    {
      MonteCarlo::LineIntersectVisit LI(tc->get<1>(),tc->get<2>());
      const Geometry::Surface* SPtr=tc->get<0>();
      
      LI.clearTrack();
      const double out=LI.getDist(SPtr);
      
      if (fabs(out-tc->get<3>())>1e-5)
	{
	  ELog::EM<<"Line :"<<tc->get<1>()<<" :: "
		  <<tc->get<2>()<<ELog::endTrace;	      
	  ELog::EM<<"Out  "<<out<<" != "<<tc->get<3>()<<ELog::endTrace;
	  ELog::EM<<"Track "<<LI.getTrack()<<ELog::endTrace;
	  ELog::EM<<"Surface :"<<*SPtr<<ELog::endTrace;	      
	  return -1;
	}
    }

  return 0;
}
  
  
