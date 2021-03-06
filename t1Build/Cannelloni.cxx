/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   t1Build/Cannelloni.cxx
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
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Triple.h"
#include "NRange.h"
#include "NList.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "surfDivide.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Sphere.h"
#include "Line.h"
#include "LineIntersectVisit.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "KGroup.h"
#include "SrcData.h"
#include "SrcItem.h"
#include "Source.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "BeamWindow.h"
#include "ProtonVoid.h"
#include "TargetBase.h"
#include "hexUnit.h"
#include "Cannelloni.h"


namespace ts1System
{

Cannelloni::Cannelloni(const std::string& Key) :
  constructSystem::TargetBase(Key,3),
  tarIndex(ModelSupport::objectRegister::Instance().cell(Key,-1,20000)),
  cellIndex(tarIndex+1),frontPlate(0),backPlate(0)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{}

Cannelloni::Cannelloni(const Cannelloni& A) : 
  constructSystem::TargetBase(A),
  tarIndex(A.tarIndex),cellIndex(A.cellIndex),
  frontPlate(A.frontPlate),backPlate(A.backPlate),
  xStep(A.xStep),yStep(A.yStep),zStep(A.zStep),
  mainLength(A.mainLength),coreRadius(A.coreRadius),
  wallThick(A.wallThick),wallClad(A.wallClad),
  frontThick(A.frontThick),voidThick(A.voidThick),
  tubeRadius(A.tubeRadius),tubeClad(A.tubeClad),
  HexHA(A.HexHA),HexHB(A.HexHB),HexHC(A.HexHC),HVec(A.HVec),
  wMat(A.wMat),taMat(A.taMat),waterMat(A.waterMat),
  targetTemp(A.targetTemp),waterTemp(A.waterTemp),
  externTemp(A.externTemp),mainCell(A.mainCell)
  /*!
    Copy constructor
    \param A :: Cannelloni to copy
  */
{}

Cannelloni&
Cannelloni::operator=(const Cannelloni& A)
  /*!
    Assignment operator
    \param A :: Cannelloni to copy
    \return *this
  */
{
  if (this!=&A)
    {
      constructSystem::TargetBase::operator=(A);
      cellIndex=A.cellIndex;
      frontPlate=A.frontPlate;
      backPlate=A.backPlate;
      xStep=A.xStep;
      yStep=A.yStep;
      zStep=A.zStep;
      mainLength=A.mainLength;
      coreRadius=A.coreRadius;
      wallThick=A.wallThick;
      wallClad=A.wallClad;
      frontThick=A.frontThick;
      voidThick=A.voidThick;
      tubeRadius=A.tubeRadius;
      tubeClad=A.tubeClad;
      HexHA=A.HexHA;
      HexHB=A.HexHB;
      HexHC=A.HexHC;
      HVec=A.HVec;
      wMat=A.wMat;
      taMat=A.taMat;
      waterMat=A.waterMat;
      targetTemp=A.targetTemp;
      waterTemp=A.waterTemp;
      externTemp=A.externTemp;
      mainCell=A.mainCell;
    }
  return *this;
}

Cannelloni*
Cannelloni::clone() const
  /*!
    Clone funciton
    \return new(this)
  */
{
  return new Cannelloni(*this);
}

Cannelloni::~Cannelloni() 
  /*!
    Destructor
  */
{}

void
Cannelloni::populate(const Simulation& System)
  /*!
    Populate all the variables
    \param System :: Simulation to use
  */
{
  const FuncDataBase& Control=System.getDataBase();

  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");

  mainLength=Control.EvalVar<double>(keyName+"MainLength");
  coreRadius=Control.EvalVar<double>(keyName+"CoreRadius");
  wallClad=Control.EvalVar<double>(keyName+"WallClad");
  wallThick=Control.EvalVar<double>(keyName+"WallThick");
  frontThick=Control.EvalVar<double>(keyName+"FrontThick");
  voidThick=Control.EvalVar<double>(keyName+"VoidThick");

  tubeRadius=Control.EvalVar<double>(keyName+"TubeRadius");
  tubeClad=Control.EvalVar<double>(keyName+"TubeClad");

  wMat=ModelSupport::EvalMat<int>(Control,keyName+"WMat");
  taMat=ModelSupport::EvalMat<int>(Control,keyName+"TaMat");
  waterMat=ModelSupport::EvalMat<int>(Control,keyName+"WaterMat");

  targetTemp=Control.EvalVar<double>(keyName+"TargetTemp");
  waterTemp=Control.EvalVar<double>(keyName+"WaterTemp");
  externTemp=Control.EvalVar<double>(keyName+"ExternTemp");


  return;
}

void
Cannelloni::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed unit for origin + xyz
  */
{
  ELog::RegMethod RegA("Cannelloni","createUnitVector");

  FixedComp::createUnitVector(FC);
  applyShift(xStep,yStep,zStep);
  // Three at 60 degrees
  HexHA=X;
  HexHB=X*cos(M_PI*60.0/180.0)+Z*sin(M_PI*60.0/180.0);
  HexHC=X*cos(M_PI*60.0/180.0)-Z*sin(M_PI*60.0/180.0);
  return;
}

void
Cannelloni::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("Cannelloni","createSurface");

  SMap.addMatch(tarIndex+1001,frontPlate);
  SMap.addMatch(tarIndex+1002,backPlate);

  // Inner Cannelloni
  ModelSupport::buildCylinder(SMap,tarIndex+7,Origin,Y,coreRadius);
  ModelSupport::buildPlane(SMap,tarIndex+1,Origin,Y);
  ModelSupport::buildPlane(SMap,tarIndex+2,Origin+Y*mainLength,Y);

  // Ta Pressure vessesl
  double T(wallClad);
  ModelSupport::buildCylinder(SMap,tarIndex+17,
			      Origin,Y,coreRadius+T);
  ModelSupport::buildPlane(SMap,tarIndex+11,Origin-Y*T,Y);
  ModelSupport::buildPlane(SMap,tarIndex+12,Origin+Y*
			   (mainLength+T),Y);
  // WAll mat
  T+=wallThick;
  ModelSupport::buildCylinder(SMap,tarIndex+27,
			      Origin,Y,coreRadius+T);
  ModelSupport::buildPlane(SMap,tarIndex+21,Origin-
			   Y*(wallClad+frontThick),Y);
  ModelSupport::buildPlane(SMap,tarIndex+22,Origin+Y*
			   (mainLength+T),Y);
  // void
  T+=voidThick;
  ModelSupport::buildCylinder(SMap,tarIndex+37,Origin,Y,
			      coreRadius+T);
  ModelSupport::buildPlane(SMap,tarIndex+31,
			   Origin-Y*(wallClad+frontThick+voidThick),Y);
  ModelSupport::buildPlane(SMap,tarIndex+32,
			   Origin+Y*(mainLength+T),Y);
  return;
}

void
Cannelloni::createObjects(Simulation& System)
  /*!
    Adds the Chip guide components
    \param System :: Simulation to create objects in
  */
{
  ELog::RegMethod RegA("Cannelloni","createObjects");
  
  // Tungsten inner core

  std::string Out;
  Out=ModelSupport::getComposite(SMap,tarIndex,"-7 1 -2");
  System.addCell(MonteCarlo::Qhull(cellIndex++,wMat,0.0,Out));
  mainCell=cellIndex-1;

  // Cladding [with front water divider]
  Out=ModelSupport::getComposite(SMap,tarIndex,"-17 11 -12 (7:-1:2) ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,taMat,0.0,Out));

  // W material
  Out=ModelSupport::getComposite(SMap,tarIndex,"-27 21 -22 (17:-11:12) ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,wMat,0.0,Out));

  // void 
  Out=ModelSupport::getComposite(SMap,tarIndex,"-37 31 -32 (27:-21:22) ");
  System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));
  
  // Set EXCLUDE:
  Out=ModelSupport::getComposite(SMap,tarIndex,"-37 31 -32");
  addOuterSurf(Out);

  addBoundarySurf(tarIndex+7);

  return;
}

void
Cannelloni::createLinks()
  /*!
    Create all the links
    \todo swap link 0 to be link 2 
  */
{
  ELog::RegMethod RegA("Cannelloni","createLinks");

  // all point out
  FixedComp::setLinkSurf(0,-SMap.realSurf(tarIndex+31));
  FixedComp::setLinkSurf(1,SMap.realSurf(tarIndex+32));
  FixedComp::setLinkSurf(2,SMap.realSurf(tarIndex+37));
		      
  FixedComp::setConnect(0,Origin-Y*(frontThick+wallClad+voidThick),-Y);
  FixedComp::setConnect(1,Origin+Y*
			(mainLength+wallClad+wallThick+voidThick),Y);
  FixedComp::setConnect(2,Origin+Z*(coreRadius+wallThick+wallClad+voidThick),Z);

  return;
}

const Geometry::Vec3D&
Cannelloni::getHexAxis(const size_t Index) const
  /*!
    Calculate the direction based on a Index 
    \param Index :: Index value [0-6]
    \return Vector for normal of plane
  */
{
  static const Geometry::Vec3D* VUnit[3]={&HexHA,&HexHC,&HexHB};
  return *(VUnit[Index % 3]);
}


void
Cannelloni::createLinkSurf()
  /*!
    Create the link surfaces
  */
{
  ELog::RegMethod RegA("Cannelloni","createLinkSurf");

  // List of surfaces constructed relative to 0,0
  std::map<int,int> surfContruct;

  std::map<int,constructSystem::hexUnit*>::iterator ac;
  std::map<int,constructSystem::hexUnit*>::iterator bc;
  int planeIndex(tarIndex+2001);
  for(ac=HVec.begin();ac!=HVec.end();ac++)
    {
      // iterate over the index set [6]
      for(size_t i=0;i<6;i++)
	{
	  constructSystem::hexUnit* APtr = ac->second;
	  if (!APtr->hasLink(i))
	    {
	      const int JA=ac->first+constructSystem::hexUnit::hexIndex(i);
	      bc=HVec.find(JA);
	      if (bc!=HVec.end())   // now construct link surface
		{
		  constructSystem::hexUnit* BPtr=bc->second; 
		  ModelSupport::buildPlane
		    (SMap,planeIndex,(APtr->getCentre()+BPtr->getCentre())/2.0,
		     getHexAxis(i));
		  int surNum=SMap.realSurf(planeIndex);
		  if (i>1 && i<5) surNum*=-1;
		  BPtr->setSurf((i+3) % 6,surNum);
		  APtr->setSurf(i,-surNum);
		  planeIndex++;
		}
	    }
	}
    }
  return;
}

void
Cannelloni::createInnerObjects(Simulation& System)
  /*!
    Build the inner objects
    \param System :: simulation
  */
{
  ELog::RegMethod RegA("Cannelloni","createInnerObject");

  int cylIndex(12000+tarIndex);
  System.removeCell(mainCell);
  const std::string endCap=
    ModelSupport::getComposite(SMap,tarIndex," 1 -2 ");

  const std::string outer=
    ModelSupport::getComposite(SMap,tarIndex," -7 ");

  std::map<int,constructSystem::hexUnit*>::const_iterator ac;
  for(ac=HVec.begin();ac!=HVec.end();ac++)
    {
      const constructSystem::hexUnit* APtr= ac->second;
      // Create Inner plane here just to help order stuff
      ModelSupport::buildCylinder(SMap,cylIndex+7,
				  APtr->getCentre(),Y,tubeRadius-tubeClad);
      ModelSupport::buildCylinder(SMap,cylIndex+8,
				  APtr->getCentre(),Y,tubeRadius-10.0*Geometry::zeroTol);
      std::string CylA=
	ModelSupport::getComposite(SMap,tarIndex,cylIndex," -7M 1 -2");
      std::string CylB=
	ModelSupport::getComposite(SMap,tarIndex,cylIndex," 7M -8M 1 -2");

      std::string Out=APtr->getInner()+
      	ModelSupport::getComposite(SMap,cylIndex," 8 ");

      if (!APtr->isComplete()) 
	{
	  CylA+=outer;
	  CylB+=outer;
	  Out+=outer;
	}
      System.addCell(MonteCarlo::Qhull(cellIndex++,wMat,0.0,CylA));
      System.addCell(MonteCarlo::Qhull(cellIndex++,taMat,0.0,CylB));
      
      System.addCell(MonteCarlo::Qhull(cellIndex++,waterMat,0.0,Out+endCap));
      cylIndex+=10;
    }
  return;
}

void 
Cannelloni::createCentres(const Geometry::Plane* PX)
  /*!
    Calculate the centres
    \param PX :: Dividing plane to test against.
  */
{
  ELog::RegMethod RegA("Cannelloni","createCentres");

  int acceptFlag(1);
  int step(0);

  HVec.clear();

  while(acceptFlag)
    {
      acceptFlag=0;
      for(int i=-step;i<=step;i++)
	for(int j=-step;j<=step;j++)
	  {
	    if (abs(i)==step || abs(j)==step)
	      {
		// Note tube *2 because separation is diameter
		const Geometry::Vec3D C=Origin+HexHA*(i*2.0*tubeRadius)+
		  HexHB*(j*2.0*tubeRadius);
		MonteCarlo::LineIntersectVisit LI(C,Y);
		const Geometry::Vec3D TPoint = LI.getPoint(PX)+
		  Y*(Geometry::zeroTol*10.0);
		const bool cFlag=isBoundaryValid(TPoint);
		if (cFlag)
		  {
		    HVec.insert(MTYPE::value_type
				((i*1000+j),
				 new constructSystem::hexUnit(i,j,C)));
		    acceptFlag=1;
		  }
	      }
	  }
      step++;
    }
  return;
}

void
Cannelloni::createInnerCells(Simulation& System)
  /*!
    Create the inner planes
    \param System :: Simulation to build
   */
{
  ELog::RegMethod RegA("Cannelloni","createInnterCells");

  const Geometry::Plane* PPtr=
    SMap.realPtr<Geometry::Plane>(tarIndex+1);

  createCentres(PPtr);
  createLinkSurf();
  createInnerObjects(System);
  return;
}

void
Cannelloni::createBeamWindow(Simulation& System)
  /*!
    Create the beamwindow if present
    \param System :: Simulation to build into
  */
{
  ELog::RegMethod RegA("Cannelloni","createBeamWindow");
  if (PLine->getVoidCell())
    {
      ModelSupport::objectRegister& OR=
	ModelSupport::objectRegister::Instance();
      
      if (!BWPtr)
	{
	  BWPtr=boost::shared_ptr<ts1System::BeamWindow>
	  (new ts1System::BeamWindow("BWindow"));
	  OR.addObject(BWPtr);
	}      
      return;
      BWPtr->addBoundarySurf(PLine->getCompContainer());
      BWPtr->setInsertCell(PLine->getVoidCell());
      BWPtr->createAll(System,*this,0);  // 0 => front face of target
    }
  return;
}


void
Cannelloni::addProtonLine(Simulation& System,
			 const attachSystem::FixedComp& refFC,
			 const long int index)
  /*!
    Add a proton void cell
    \param System :: Simualation
    \param refFC :: reflector edge
    \param index :: Index of the proton cutting surface [6 typically (-7)]
   */
{
  ELog::RegMethod RegA("t1PlateTarget","addProtonLine");

  // 0 ::  front fact of target
  PLine->createAll(System,*this,0,refFC,index);
  createBeamWindow(System);
  System.populateCells();
  System.createObjSurfMap();
  return;
}

  
void
Cannelloni::createAll(Simulation& System,
		       const attachSystem::FixedComp& FC)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: Fixed Component for origin
  */
{
  ELog::RegMethod RegA("Cannelloni","createAll");

  populate(System);
  createUnitVector(FC);
  createSurfaces();
  createObjects(System);
  createLinks();
  createBeamWindow(System);
  createInnerCells(System);
  insertObjects(System);


  return;
}


  
}  // NAMESPACE TMRsystem
