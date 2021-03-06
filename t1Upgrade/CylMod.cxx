/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   t1Upgrade/CylMod.cxx
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
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "Quadratic.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "support.h"
#include "stringCombine.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "LayerComp.h"
#include "ConicInfo.h"
#include "ModBase.h"
#include "CylMod.h"

namespace ts1System
{

CylMod::CylMod(const std::string& Key) :
  ModBase(Key,4)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}

CylMod::CylMod(const CylMod& A) : 
  ModBase(A),
  nLayers(A.nLayers),radius(A.radius),
  height(A.height),mat(A.mat)
  /*!
    Copy constructor
    \param A :: CylMod to copy
  */
{}

CylMod&
CylMod::operator=(const CylMod& A)
  /*!
    Assignment operator
    \param A :: CylMod to copy
    \return *this
  */
{
  if (this!=&A)
    {
      ModBase::operator=(A),
      cellIndex=A.cellIndex;
      nLayers=A.nLayers;
      radius=A.radius;
      height=A.height;
      mat=A.mat;
    }
  return *this;
}

CylMod::~CylMod()
  /*!
    Destructor
  */
{}

CylMod*
CylMod::clone() const
  /*!
    Clone copy constructor
    \return copy of this
  */
{
  return new CylMod(*this);
}

void
CylMod::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("CylMod","populate");

    // Master values
  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xyAngle=Control.EvalVar<double>(keyName+"XYangle");
  zAngle=Control.EvalVar<double>(keyName+"Zangle");

  double R,H,T;
  int M;
  nLayers=Control.EvalVar<size_t>(keyName+"NLayers");   
  for(size_t i=0;i<=nLayers;i++)
    {
      if (i)
	{
	  H+=2.0*Control.EvalVar<double>
	    (StrFunc::makeString(keyName+"HGap",i));   
	  R+=Control.EvalVar<double>
	    (StrFunc::makeString(keyName+"RadGap",i));   
	  M=ModelSupport::EvalMat<int>
	    (Control,StrFunc::makeString(keyName+"Material",i));   
	  T=(!M) ? 0.0 : 
	    Control.EvalVar<double>(StrFunc::makeString(keyName+"Temp",i)); 
	}
      else
	{
	  H=Control.EvalVar<double>(keyName+"Height");   
	  R=Control.EvalVar<double>(keyName+"Radius");   
	  M=ModelSupport::EvalMat<int>(Control,keyName+"Mat");   
	  T=Control.EvalVar<double>(keyName+"Temp");   
	}
      radius.push_back(R);
      height.push_back(H);
      mat.push_back(M);
      temp.push_back(T);
    }
  nConic=Control.EvalVar<size_t>(keyName+"NConic");
  for(size_t i=0;i<nConic;i++)
    {
      const std::string KN=keyName+StrFunc::makeString("Conic",i+1);
      const Geometry::Vec3D C=
	Control.EvalVar<Geometry::Vec3D>(KN+"Cent");
      const Geometry::Vec3D A=
	Control.EvalVar<Geometry::Vec3D>(KN+"Axis");
      const double W=
	Control.EvalVar<double>(KN+"Wall");
      const int WM=
	Control.EvalVar<int>(KN+"WallMat");
      const double ang=
	Control.EvalVar<double>(KN+"Angle");
      Conics.push_back(ConicInfo(C,A,ang,W,WM));
    }
  
  return;
}

void
CylMod::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("CylMod","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);

  applyShift(xStep,yStep,zStep);
  applyAngleRotate(xyAngle,zAngle);


  return;
}

void
CylMod::createSurfaces()
  /*!
    Create planes for the silicon and Polyethene layers
  */
{
  ELog::RegMethod RegA("CylMod","createSurfaces");

  // Divide plane
  ModelSupport::buildPlane(SMap,modIndex+1,Origin,X);  
  ModelSupport::buildPlane(SMap,modIndex+2,Origin,Y);  

  int SI(modIndex);
  for(size_t i=0;i<=nLayers;i++)
    {
      ModelSupport::buildCylinder(SMap,SI+7,Origin,Z,radius[i]);  
      ModelSupport::buildPlane(SMap,SI+5,Origin-Z*height[i]/2.0,Z);  
      ModelSupport::buildPlane(SMap,SI+6,Origin+Z*height[i]/2.0,Z);  
      SI+=10;
    }
  // CONICS
  SI=modIndex+500;
  for(size_t i=0;i<nConic;i++)
    {
      const Geometry::Vec3D Pt=Conics[i].getCent(X,Y,Z);
      const Geometry::Vec3D CY=Conics[i].getAxis(X,Y,Z);
      const int cylFlag=Conics[i].isCyl();
      const double rangle=Conics[i].getAngle();
      ModelSupport::buildPlane(SMap,SI+1,Origin+Pt,CY);
      if (!cylFlag)
	{
	  const double wallThick=Conics[i].getWall()/cos(M_PI*rangle/180.0);
	  // Note dividing surface
	  ModelSupport::buildCone(SMap,SI+7,Origin+Pt,CY,rangle);
	  if (Conics[i].getWall()>Geometry::zeroTol)
	    {
	      ModelSupport::buildCone(SMap,SI+17,
				      Origin+Pt-CY*wallThick,CY,rangle);
	      ModelSupport::buildPlane(SMap,SI+11,Origin+Pt-CY*wallThick,CY);
	    }
	}
      else
	{
	  const double wallThick=Conics[i].getWall();
	  ModelSupport::buildCylinder(SMap,SI+7,Origin+Pt,CY,rangle);
	  ModelSupport::buildPlane(SMap,SI+1,Origin+Pt,CY);
	  if (Conics[i].getWall()>Geometry::zeroTol)
	    {
	      ModelSupport::buildCylinder(SMap,SI+17,Origin+Pt,CY,
					  rangle+wallThick);
	      ModelSupport::buildPlane(SMap,SI+11,Origin+Pt-CY*wallThick,CY);
	    }
	}

      SI+=100;
    }

  return; 
}

void
CylMod::createObjects(Simulation& System)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
  */
{
  ELog::RegMethod RegA("CylMod","createObjects");
  const FuncDataBase& Control =System.getDataBase();
  const int cylFlag=Control.EvalDefVar<int>(keyName+"CylFlag",0);
  ELog::EM<<"cylFlag = "<<cylFlag<<ELog::endDebug;
  std::string Out;
  // First make conics:
  int CI(modIndex+500);
  HeadRule OutUnit;
  for(size_t i=0;i<nConic;i++)
    {
      Out=ModelSupport::getComposite(SMap,modIndex,CI," -7 5 -6 -7M 1M");
      System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));
      if (Conics[i].getWall()>Geometry::zeroTol)
	{
	  Out=ModelSupport::getComposite(SMap,modIndex,CI,
					 " -7 5 -6 (7M:-1M) -17M 11M");
	  System.addCell(MonteCarlo::Qhull(cellIndex++,
					   Conics[i].getMat(),
					   temp[0],Out));
	  Out=ModelSupport::getComposite(SMap,modIndex,CI," -7 5 -6 -17M 11M ");
	}
      OutUnit.addUnion(Out);
      CI+=100;
    }

  OutUnit.makeComplement();

  int SI(modIndex);
  for(size_t i=0;i<=nLayers;i++)
    {
      Out=ModelSupport::getComposite(SMap,SI," -7 5 -6 ");

      if (i==nLayers) addOuterSurf(Out);
      if (i)
	Out+=ModelSupport::getComposite(SMap,SI-10," (7:-5:6) ");
      else
	Out+=OutUnit.display();
      System.addCell(MonteCarlo::Qhull(cellIndex++,mat[i],temp[i],Out));
      SI+=10;
    }
  return; 
}

void
CylMod::createLinks()
  /*!
    Creates a full attachment set
  */
{  
  ELog::RegMethod RegA("CylMod","createLinks");

  const int SI(modIndex+static_cast<int>(nLayers)*10);
  
  FixedComp::setConnect(0,Origin-Y*radius[nLayers],-Y);
  FixedComp::setLinkSurf(0,SMap.realSurf(SI+7));
  FixedComp::addLinkSurf(0,-SMap.realSurf(modIndex+2));
  
  FixedComp::setConnect(1,Origin+Y*radius[nLayers],Y);
  FixedComp::setLinkSurf(1,SMap.realSurf(SI+7));
  FixedComp::addLinkSurf(1,SMap.realSurf(modIndex+2));
  
  FixedComp::setConnect(2,Origin-Z*(height[nLayers]/2.0),-Z);
  FixedComp::setLinkSurf(2,-SMap.realSurf(SI+5));
  
  FixedComp::setConnect(3,Origin+Z*(height[nLayers]/2.0),Z);
  FixedComp::setLinkSurf(3,SMap.realSurf(SI+6));

  return;
}


Geometry::Vec3D
CylMod::getSurfacePoint(const size_t layerIndex,
			const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link point
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-6]
    \return Surface point
  */
{
  ELog::RegMethod RegA("CylMod","getSurfacePoint");

  if (sideIndex>5) 
    throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

  // Modification map:
  switch(sideIndex)
    {
    case 0:
      return Origin-Y*radius[layerIndex];
    case 1:
      return Origin+Y*radius[layerIndex];
    case 2:
      return Origin-X*radius[layerIndex];
    case 3:
      return Origin+X*radius[layerIndex];
    case 4:
      return Origin-Z*(height[layerIndex]/2.0);
    case 5:
      return Origin+Z*(height[layerIndex]/2.0);
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

std::string
CylMod::getLayerString(const size_t sideIndex,
		       const size_t layerIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
{
  ELog::RegMethod RegA("CylMod","getLayerString");

  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layer");

  const int SI(modIndex+static_cast<int>(layerIndex)*10);
  std::ostringstream cx;
  switch(sideIndex)
    {
    case 0:
    case 1:
    case 2:
    case 3:
      cx<<" "<<SMap.realSurf(SI+7)<<" ";
      return cx.str();
    case 4:
      cx<<" "<<-SMap.realSurf(SI+5)<<" ";
      return cx.str();
    case 5:
      cx<<" "<<SMap.realSurf(SI+6)<<" ";
      return cx.str();
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

int
CylMod::getLayerSurf(const size_t layerIndex,
		     const size_t sideIndex) const
  /*!
    Given a side and a layer calculate the link surf
    \param sideIndex :: Side [0-5]
    \param layerIndex :: layer, 0 is inner moderator [0-4]
    \return Surface string
  */
{
  ELog::RegMethod RegA("CylMod","getLayerSurf");

  if (layerIndex>=nLayers) 
    throw ColErr::IndexError<size_t>(layerIndex,nLayers,"layerIndex");
  
  const int SI(modIndex+static_cast<int>(layerIndex)*10);
  switch(sideIndex)
    {
    case 0:
    case 1:
    case 2:
    case 3:
      return SMap.realSurf(SI+7);
    case 4:
      return -SMap.realSurf(SI+5);
    case 5:
      return SMap.realSurf(SI+6);
    }
  throw ColErr::IndexError<size_t>(sideIndex,5,"sideIndex ");
}

void
CylMod::createAll(Simulation& System,
		     const attachSystem::FixedComp& FC)
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
   */
{
  ELog::RegMethod RegA("CylMod","createAll");
  populate(System.getDataBase());

  createUnitVector(FC);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
