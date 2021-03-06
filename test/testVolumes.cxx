/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   test/testVolumes.cxx
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
#include <cmath>
#include <complex> 
#include <vector>
#include <list> 
#include <map> 
#include <set>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <boost/functional.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>

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
#include "version.h"
#include "Element.h"
#include "MapSupport.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Triple.h"
#include "NList.h"
#include "NRange.h"
#include "Tally.h"
#include "cellFluxTally.h"
#include "pointTally.h"
#include "heatTally.h"
#include "tallyFactory.h"
#include "Transform.h"
#include "Surface.h"
#include "surfIndex.h"
#include "Quadratic.h"
#include "surfaceFactory.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FItem.h"
#include "FuncDataBase.h"
#include "SurInter.h"
#include "BnId.h"
#include "Acomp.h"
#include "Algebra.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "RemoveCell.h"
#include "WForm.h"
#include "weightManager.h"
#include "ObjSurfMap.h"
#include "ObjTrackItem.h"
#include "SrcData.h"
#include "SrcItem.h"
#include "Source.h"
#include "ReadFunctions.h"
#include "surfRegister.h"
#include "ModelSupport.h"
#include "neutron.h"
#include "inputParam.h"
#include "Simulation.h"
#include "volUnit.h"
#include "VolSum.h"
#include "Volumes.h"

#include "testFunc.h"
#include "testVolumes.h"

using namespace ModelSupport;

testVolumes::testVolumes() 
  /*!
    Constructor
  */
{
  initSim();
}

testVolumes::~testVolumes() 
  /*!
    Destructor
  */
{}

void
testVolumes::initSim()
  /*!
    Set all the objects in the simulation:
  */
{
  ASim.resetAll();
  createSurfaces();
  createObjects();
  ASim.createObjSurfMap();
  return;
}

void 
testVolumes::createSurfaces()
  /*!
    Create the surface list
   */
{
  ELog::RegMethod RegA("testVolumes","createSurfaces");

  ModelSupport::surfIndex& SurI=ModelSupport::surfIndex::Instance();
  
  // First box :
  SurI.createSurface(1,"px -2");
  SurI.createSurface(2,"px 2");
  SurI.createSurface(3,"py -1");
  SurI.createSurface(4,"py 1");
  SurI.createSurface(5,"pz -1");
  SurI.createSurface(6,"pz 1");

  // Second box :
  SurI.createSurface(11,"px -3");
  SurI.createSurface(12,"px 3");
  SurI.createSurface(13,"py -3");
  SurI.createSurface(14,"py 3");
  SurI.createSurface(15,"pz -3");
  SurI.createSurface(16,"pz 3");

  // Far box :
  SurI.createSurface(21,"px 10");
  SurI.createSurface(22,"px 15");

  // Sphere :
  SurI.createSurface(100,"so 25");
  // Sphere :
  SurI.createSurface(101,"so 6.0");
  SurI.createSurface(102,"s 7.0 0.0 0.0 3.0");
  
  return;
}
  
void
testVolumes::createObjects()
  /*!
    Create Object for test
   */
{
  std::string Out;
  int cellIndex(1);
  const int surIndex(0);

  Out=ModelSupport::getComposite(surIndex,"100");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));      // Outside void Void

  // Out=ModelSupport::getComposite(surIndex,"-102");
  // ASim.addCell(MonteCarlo::Qhull(cellIndex++,5,0.0,Out));   
  // Out=ModelSupport::getComposite(surIndex,"-101");
  // ASim.addCell(MonteCarlo::Qhull(cellIndex++,5,0.0,Out));   
  // Out=ModelSupport::getComposite(surIndex,"-100 101 102");
  // ASim.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));   

  //  Out=ModelSupport::getComposite(surIndex,"1 -2 3 -4 5 -6");
  //  ASim.addCell(MonteCarlo::Qhull(cellIndex++,3,0.0,Out));      // steel object

  Out=ModelSupport::getComposite(surIndex,"-101");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,78,0.0,Out));      // steel object

  // Out=ModelSupport::getComposite(surIndex,"-100 (-1:2:-3:4:-5:6)");
  Out=ModelSupport::getComposite(surIndex,"-100 101");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));      // steel object

  ASim.populateCells();
  return;

  Out=ModelSupport::getComposite(surIndex,"11 -12 13 -14 15 -16"
				 " (-1:2:-3:4:-5:6) ");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,5,0.0,Out));      // Al container

  Out=ModelSupport::getComposite(surIndex,"21 -22 3 -4 5 -6");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,8,0.0,Out));      // Gd box 

  Out=ModelSupport::getComposite(surIndex,"-100 (-11:12:-13:14:-15:16)"
				 " #4");
  ASim.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));      // Void
  
  ASim.removeComplements();

  return;
}
  

int 
testVolumes::applyTest(const int extra)
  /*!
    Applies all the tests and returns 
    the error number
    \param extra :: Test number to run
    \retval -1 : SetObject 
    \retval 0 : All succeeded
  */
{
  ELog::RegMethod RegA("testVolumes","applyTest");
  TestFunc::regSector("testVolumes");

  typedef int (testVolumes::*testPtr)();
  testPtr TPtr[]=
    {
      &testVolumes::testPointVolume,
      &testVolumes::testVolume
    };
  const std::string TestName[]=
    {
      "PointVolume",
      "Volume"
    };
  
  const int TSize(sizeof(TPtr)/sizeof(testPtr));
  if (!extra)
    {
      std::ios::fmtflags flagIO=std::cout.setf(std::ios::left);
      for(int i=0;i<TSize;i++)
	std::cout<<std::setw(30)<<TestName[i]<<"("<<i+1<<")"<<std::endl;
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
testVolumes::testVolume()
  /*!
    Tracks a neutron through the system
    \return 0 on success and -1 on error
  */
{
  ELog::RegMethod RegA("testVolumes","testVolume");
  double V;


  VolSum VTallyX(Geometry::Vec3D(0,0,-4),8.0);
  //  VTally.populate(ASim);
  VTallyX.addTallyCell(4,2);
  //  VTally.addTallyCell(14,3);
  VTallyX.run(ASim,80000);
  V=VTallyX.calcVolume(4);
  ELog::EM<<"CalcValue(4X) == "<<V<<ELog::endTrace;
  ELog::EM<<" ========= "<<V<<ELog::endTrace<<
    ELog::endTrace;
  VolSum VTally(Geometry::Vec3D(0,0,0),8.0);
  //  VTally.populate(ASim);
  VTally.addTallyCell(4,2);
  //  VTally.addTallyCell(14,3);
  VTally.run(ASim,80000);
  V=VTally.calcVolume(4);
  ELog::EM<<"CalcValue(4) == "<<V<<ELog::endTrace;
  ELog::EM<<" ========= "<<V<<ELog::endTrace<<
    ELog::endTrace;




  return 0;
  VTally.run(ASim,20000);
  ELog::EM<<"CalcValue(4) == "<<VTally.calcVolume(4)<<ELog::endTrace;
  VTally.run(ASim,200000);
  ELog::EM<<"CalcValue(4) == " <<VTally.calcVolume(4)<<ELog::endTrace;
  //  VTally.run(ASim,2000000);
  //  ELog::EM<<"CalcValue(4) == "<<VTally.calcVolume(4)<<ELog::endTrace;

  //  ELog::EM<<"CalcValue(14) == "<<VTally.calcVolume(14)<<ELog::endTrace;
  
  return 0;
}

int
testVolumes::testPointVolume()
  /*!
    Tracks a neutron through the system
    \return 0 on success and -1 on error
  */
{
  ELog::RegMethod RegA("testVolumes","testVolume");
  double V;


  VolSum VTallyX(Geometry::Vec3D(0,0,0),8.0);
  //  VTally.populate(ASim);
  VTallyX.addTallyCell(4,2);
  VTallyX.addTallyCell(5,3);
  //  VTally.addTallyCell(14,3);
  VTallyX.pointRun(ASim,800000);
  V=VTallyX.calcVolume(4);
  double Vx=VTallyX.calcVolume(5);
  ELog::EM<<"CalcValue(4X) == "<<V<<ELog::endTrace;
  ELog::EM<<" ========= "<<V<<ELog::endTrace<<
    ELog::endTrace;




  
  return 0;
}

