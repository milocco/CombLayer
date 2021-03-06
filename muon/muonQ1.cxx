/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   muon/muonQ1.cxx
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
#include <iterator>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>

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
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Rules.h"
#include "Convex.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "SimProcess.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "ContainedComp.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "World.h"
#include "muonQ1.h"

namespace muSystem
{

muonQ1::muonQ1(const std::string& Key)  : 
  attachSystem::FixedComp(Key,6),attachSystem::ContainedComp(),
  muQ1Index(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(muQ1Index+1)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Key to use
  */
{}


muonQ1::~muonQ1() 
  /*!
    Destructor
  */
{}

void
muonQ1::populate(const Simulation& System)
  /*!
    Populate all the variables
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("muonQ1","populate");

  const FuncDataBase& Control=System.getDataBase();

  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xAngle=Control.EvalVar<double>(keyName+"Xangle");
  yAngle=Control.EvalVar<double>(keyName+"Yangle");  
  zAngle=Control.EvalVar<double>(keyName+"Zangle");

  xSize=Control.EvalVar<double>(keyName+"XSize");
  ySize=Control.EvalVar<double>(keyName+"YSize");
  zSize=Control.EvalVar<double>(keyName+"ZSize");
  cutLenOut=Control.EvalVar<double>(keyName+"CutLenOut");
  cutLenMid=Control.EvalVar<double>(keyName+"CutLenMid");

  steelThick=Control.EvalVar<double>(keyName+"SteelThick");
  copperThick=Control.EvalVar<double>(keyName+"CopperThick");
  copperYSize=Control.EvalVar<double>(keyName+"CopperYSize");    

  insertSize=Control.EvalVar<double>(keyName+"InsertSize");    
  insertThick=Control.EvalVar<double>(keyName+"InsertThick");    
      
  steelMat=ModelSupport::EvalMat<int>(Control,keyName+"SteelMat");    
  copperMat=ModelSupport::EvalMat<int>(Control,keyName+"CopperMat");
  insertMat=ModelSupport::EvalMat<int>(Control,keyName+"InsertMat");
         
  return;
}

void
muonQ1::createUnitVector()
  /*!
    Create the unit vectors
  */
{
  ELog::RegMethod RegA("muonQ1","createUnitVector");

  attachSystem::FixedComp::createUnitVector(World::masterOrigin());
  applyShift(xStep,yStep,zStep);
  applyAngleRotate(0,0,zAngle);  
  applyAngleRotate(0,yAngle,0);  
//  applyShift(xStep,yStep,zStep);


  
  return;
}

void
muonQ1::createSurfaces()
  /*!
    Create all the surfaces
  */
{
  ELog::RegMethod RegA("muonQ1","createSurface");

  // Outer layer:
  ModelSupport::buildPlane(SMap,muQ1Index+1,Origin-Y*ySize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+2,Origin+Y*ySize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+3,Origin-X*xSize/2.0,X);
  ModelSupport::buildPlane(SMap,muQ1Index+4,Origin+X*xSize/2.0,X);
  ModelSupport::buildPlane(SMap,muQ1Index+5,Origin-Z*zSize/2.0,Z);
  ModelSupport::buildPlane(SMap,muQ1Index+6,Origin+Z*zSize/2.0,Z);
         // Corners:
  ModelSupport::buildPlane(SMap,muQ1Index+101,Origin
			   -X*xSize/2.0-Z*(zSize/2.0-cutLenOut),-X-Z);
  ModelSupport::buildPlane(SMap,muQ1Index+102,Origin
			   -X*xSize/2.0+Z*(zSize/2.0-cutLenOut),-X+Z);
  ModelSupport::buildPlane(SMap,muQ1Index+103,Origin
			   +X*xSize/2.0-Z*(zSize/2.0-cutLenOut),X-Z);
  ModelSupport::buildPlane(SMap,muQ1Index+104,Origin
			   +X*xSize/2.0+Z*(zSize/2.0-cutLenOut),X+Z);

  // Steel layer:
  ModelSupport::buildPlane(SMap,muQ1Index+11,Origin-Y*copperYSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+12,Origin+Y*copperYSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+13,Origin-X*(xSize/2.0-steelThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+14,Origin+X*(xSize/2.0-steelThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+15,Origin-Z*(zSize/2.0-steelThick),Z);
  ModelSupport::buildPlane(SMap,muQ1Index+16,Origin+Z*(zSize/2.0-steelThick),Z);
         // Corners:
  ModelSupport::buildPlane(SMap,muQ1Index+111,Origin
			   -X*(xSize/2.0-steelThick)-Z*(zSize/2.0-steelThick-cutLenMid),-X-Z);
  ModelSupport::buildPlane(SMap,muQ1Index+112,Origin
			   -X*(xSize/2.0-steelThick)+Z*(zSize/2.0-steelThick-cutLenMid),-X+Z);
  ModelSupport::buildPlane(SMap,muQ1Index+113,Origin
			   +X*(xSize/2.0-steelThick)-Z*(zSize/2.0-steelThick-cutLenMid),X-Z);
  ModelSupport::buildPlane(SMap,muQ1Index+114,Origin
			   +X*(xSize/2.0-steelThick)+Z*(zSize/2.0-steelThick-cutLenMid),X+Z);

  // Coper layer:
//  ModelSupport::buildPlane(SMap,muQ1Index+21,Origin-Y*copperYSize/2.0,Y);
//  ModelSupport::buildPlane(SMap,muQ1Index+22,Origin+Y*copperYSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+23,Origin-X*(xSize/2.0-steelThick-copperThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+24,Origin+X*(xSize/2.0-steelThick-copperThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+25,Origin-Z*(zSize/2.0-steelThick-copperThick),Z);
  ModelSupport::buildPlane(SMap,muQ1Index+26,Origin+Z*(zSize/2.0-steelThick-copperThick),Z);

  // Inserts - top & bottom:
  ModelSupport::buildPlane(SMap,muQ1Index+31,Origin-Y*insertSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+32,Origin+Y*insertSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+33,Origin-X*insertSize/2.0,X);
  ModelSupport::buildPlane(SMap,muQ1Index+34,Origin+X*insertSize/2.0,X);
  ModelSupport::buildPlane(SMap,muQ1Index+35,Origin-Z*(zSize/2.0-steelThick-copperThick-insertThick),Z);
  ModelSupport::buildPlane(SMap,muQ1Index+36,Origin+Z*(zSize/2.0-steelThick-copperThick-insertThick),Z);
  
   // Inserts - left & right:
//  ModelSupport::buildPlane(SMap,muQ1Index+31,Origin-Y*insertSize/2.0,Y);
//  ModelSupport::buildPlane(SMap,muQ1Index+32,Origin+Y*insertSize/2.0,Y);
  ModelSupport::buildPlane(SMap,muQ1Index+43,Origin-X*(xSize/2.0-steelThick-copperThick-insertThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+44,Origin+X*(xSize/2.0-steelThick-copperThick-insertThick),X);
  ModelSupport::buildPlane(SMap,muQ1Index+45,Origin-Z*insertSize/2.0,Z);
  ModelSupport::buildPlane(SMap,muQ1Index+46,Origin+Z*insertSize/2.0,Z); 
  
  return;
}

void
muonQ1::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=muQ1Index+1;i<cellIndex;i++)
    CC.addInsertCell(i);
  return;
}

void
muonQ1::createObjects(Simulation& System)
  /*!
    Adds the Chip guide components
    \param System :: Simulation to create objects in
   */
{
  ELog::RegMethod RegA("muonQ1","createObjects");
  
  std::string Out;
  std::string Out1;
  std::string Out2;
  std::string Out3;
  std::string Out4;
  std::string Out5;
    
    // Steel
  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "1 -2 3 -4 5 -6 -101 -102 -103 -104 ");
  addOuterSurf(Out);
  addBoundarySurf(Out);
  Out1=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-13:14:-15:16:111:112:113:114)");  
  System.addCell(MonteCarlo::Qhull(cellIndex++,steelMat,0.0,Out+Out1));  

    // Copper
  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "11 -12 13 -14 15 -16 -111 -112 -113 -114 ");
  addOuterUnionSurf(Out);
  addBoundaryUnionSurf(Out);
  Out1=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-23:24:-25:26)");
  Out2=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-31:32:-33:34:-36:16)");  	
  Out3=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-31:32:-33:34:35:-15)");
  Out4=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-31:32:-13:43:-45:46)"); 
  Out5=ModelSupport::getComposite(SMap,muQ1Index,
				 "(-31:32:-44:14:-45:46)");  					  					   				   
  System.addCell(MonteCarlo::Qhull(cellIndex++,copperMat,0.0,Out+Out1+Out2+Out3+Out4+Out5));  


    // Inserts
  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "31 -32 33 -34 36 -16 ");    
  System.addCell(MonteCarlo::Qhull(cellIndex++,insertMat,0.0,Out));    

  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "31 -32 33 -34 -35 15 ");    
  System.addCell(MonteCarlo::Qhull(cellIndex++,insertMat,0.0,Out));    

  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "31 -32 13 -43 45 -46 ");    
  System.addCell(MonteCarlo::Qhull(cellIndex++,insertMat,0.0,Out));    

  Out=ModelSupport::getComposite(SMap,muQ1Index,
				 "31 -32 44 -14 45 -46 ");    
  System.addCell(MonteCarlo::Qhull(cellIndex++,insertMat,0.0,Out));    

    // Void
  Out=ModelSupport::getComposite(SMap,muQ1Index,"11 -12 23 -24 25 -26 ");
  
  System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out+Out2+Out3+Out4+Out5));    
  return;
}


void
muonQ1::createLinks()
  /*!
    Create links
   */
{
  ELog::RegMethod RegA("muonQ1","createLinks");

  FixedComp::setLinkSurf(0,-SMap.realSurf(muQ1Index+1));
  FixedComp::setLinkSurf(1,SMap.realSurf(muQ1Index+2));
  FixedComp::setLinkSurf(2,-SMap.realSurf(muQ1Index+3));
  FixedComp::setLinkSurf(3,SMap.realSurf(muQ1Index+4));
  FixedComp::setLinkSurf(4,-SMap.realSurf(muQ1Index+5));
  FixedComp::setLinkSurf(5,SMap.realSurf(muQ1Index+6));

  return;
}

void
muonQ1::createAll(Simulation& System)

  /*!
    Global creation of the hutch
    \param System :: Simulation to add vessel to
    \param FC :: Fixed Component to place object within
  */
{
  ELog::RegMethod RegA("muonQ1","createAll");
  populate(System);
  createUnitVector();
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);
  return;
}
  
}  // NAMESPACE shutterSystem
