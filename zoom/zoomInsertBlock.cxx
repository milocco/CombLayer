/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   zoom/zoomInsertBlock.cxx
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
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <boost/format.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Triple.h"
#include "NRange.h"
#include "NList.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "KGroup.h"
#include "Source.h"
#include "shutterBlock.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinearComp.h"
#include "ContainedComp.h"
#include "zoomInsertBase.h"
#include "zoomInsertBlock.h"

namespace shutterSystem
{


zoomInsertBlock::zoomInsertBlock(const int N,const int SN,
				 const std::string& Key) :
  zoomInsertBase(N,SN,Key)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param N :: Index value of block
    \param SN :: surface Index value
    \param Key :: Name for item in search
  */
{}

zoomInsertBlock::zoomInsertBlock(const zoomInsertBlock& A) : 
  zoomInsertBase(A),
  hGap(A.hGap),vGap(A.vGap)
  /*!
    Copy constructor
    \param A :: zoomInsertBlock to copy
  */
{}

zoomInsertBlock&
zoomInsertBlock::operator=(const zoomInsertBlock& A)
  /*!
    Assignment operator
    \param A :: zoomInsertBlock to copy
    \return *this
  */
{
  if (this!=&A)
    {
      zoomInsertBase::operator=(A);
      hGap=A.hGap;
      vGap=A.vGap;
    }
  return *this;
}

zoomInsertBlock*
zoomInsertBlock::clone() const
  /*!
    Clone constructor
    \return copy of this
  */
{
  return new zoomInsertBlock(*this);
}

void
zoomInsertBlock::populate(const Simulation& System,
			  const zoomInsertBase* sndBase)
  /*!
    Populate all the variables
    Requires that unset values are copied from previous block
    \param System :: Simulation to use
    \param sndBlock :: Standard block to use
  */
{
  ELog::RegMethod RegA("zoomInsertBlock","populate");
  const FuncDataBase& Control=System.getDataBase();

  const int Size(9);
  const int commonSize(7);
  const char* sndKey[Size]=
    {"FStep","CentX","CentZ","Len","Width","Height","Mat","HGap","VGap"};
  
  const zoomInsertBlock* sndBlock=
    dynamic_cast<const zoomInsertBlock*>(sndBase);

  for(int i=0;i<Size;i++)
    {
      std::ostringstream cx;
      cx<<keyName<<blockIndex+1<<sndKey[i];
      if (Control.hasVariable(cx.str()))
	setVar(Control,i,cx.str());
      else if (sndBlock)
	setVar(i,sndBlock->getVar(i));	
      else if (sndBase && i<=commonSize)
	 setVar(i,sndBase->getVar(i));	
      else 
	{
	  ELog::EM<<"sndBase == "<<sndBase->typeName()<<ELog::endCrit;
	  ELog::EM<<"Failed to connect on first component:"
		  <<blockIndex+1<<" :: "<<cx.str()<<ELog::endErr;
	}
    }
  populated|=1;
  return;
}

void
zoomInsertBlock::setVar(const int Item,const double V)
  /*!
    Given a value set the item
    \param Item :: Index value to variable
    \param V :: Value
  */
{
  switch(Item)
    {
    case 0:
      fStep=V;
      return;
    case 1:
      centX=V;
      return;
    case 2:
      centZ=V;
      return;
    case 3:
      length=V;
      return;
    case 4:
      width=V;
      return;
    case 5:
      height=V;
      return;
    case 6:
      matN=static_cast<int>(V);
      return;
    case 7:
      hGap=V;
      return;
    case 8:
      vGap=V;
      return;
    default:
      throw ColErr::IndexError<int>(Item,8,"Item");
    }
  return;
}

void
zoomInsertBlock::setVar(const FuncDataBase& Control,
			const int Item,const std::string& VarStr)
  /*!
    Convert a string into an item value
    \param Control :: DataBase to get value from
    \param Item :: Index value to variable
    \param VarStr :: String to convert into values
  */
{
  ELog::RegMethod RegA("zoomInsertBlock","setVar<Control>");
  setVar(Item,Control.EvalVar<double>(VarStr));
  return;
}

double
zoomInsertBlock::getVar(const int Item) const
  /*!
    Get the value based on an index reference
    \param Item :: Index value to variable
    \return value
  */
{
  ELog::RegMethod RegA("zoomInsertBlock","getVar");
  
  switch(Item)
    {
    case 0:
      return fStep;
    case 1:
      return centX;
    case 2:
      return centZ;
    case 3:
      return length;
    case 4:
      return width;
    case 5:
      return height;
    case 6:
      return matN;
    case 7:
      return hGap;
    case 8:
      return vGap;
    }
  throw ColErr::IndexError<int>(Item,8,"Item");
}

void
zoomInsertBlock::createSurfaces(const int startSurf)
  /*!
    Creates/duplicates the surfaces for this block
    \param startSurf :: First surface to use
  */
{  
  ELog::RegMethod RegA("zoomInsertBlock","createSurface");

  if (startSurf)
    SMap.addMatch(surfIndex+1,startSurf);
  else
    ModelSupport::buildPlane(SMap,surfIndex+1,Origin,Y);
  
  // Outer Surface 
  ModelSupport::buildPlane(SMap,surfIndex+2,Origin+Y*length,Y);
  ModelSupport::buildPlane(SMap,surfIndex+3,Origin-X*width,X);
  ModelSupport::buildPlane(SMap,surfIndex+4,Origin+X*width,X);
  ModelSupport::buildPlane(SMap,surfIndex+5,Origin-Z*height,Z);
  ModelSupport::buildPlane(SMap,surfIndex+6,Origin+Z*height,Z);
  
  // Inner surface
  ModelSupport::buildPlane(SMap,surfIndex+13,
			   beamOrigin-beamX*(hGap/2.0-centX),beamX);

  ModelSupport::buildPlane(SMap,surfIndex+14,
			   beamOrigin+beamX*(hGap/2.0+centX),beamX);

  ModelSupport::buildPlane(SMap,surfIndex+15,
			   beamOrigin-beamZ*(vGap/2.0-centZ),beamZ);
  ModelSupport::buildPlane(SMap,surfIndex+16,
			   beamOrigin+beamZ*(vGap/2.0+centZ),beamZ);
  
  return;
}

void
zoomInsertBlock::createObjects(Simulation& System,
			       const std::string& fSurf,
			       const std::string& bSurf)
  /*!
    Create the objects
    \param System :: Simulation
    \param fSurf :: front surface 
    \param bSurf :: back Surface
  */
{
  ELog::RegMethod RegA("zoomInsertBlock","createObjects");

  std::string frontBack=fSurf.empty() ? 
    ModelSupport::getComposite(SMap,surfIndex,"1 ") : fSurf;
  frontBack+=bSurf.empty() ? 
    ModelSupport::getComposite(SMap,surfIndex,"-2 ") : bSurf;

  std::string Out;
  Out=ModelSupport::getComposite(SMap,surfIndex,"3 -4 5 -6 ");
  addOuterSurf(Out);

  // Centre void
  Out=ModelSupport::getComposite(SMap,surfIndex,"13 -14 15 -16 ")+
    frontBack;
  System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));

  // Outer metal
  Out=ModelSupport::getComposite(SMap,surfIndex,
				 "3 -4 5 -6 (-13:14:-15:16) ")+
    frontBack;
  System.addCell(MonteCarlo::Qhull(cellIndex++,matN,0.0,Out));
  
  return;
}

int
zoomInsertBlock::exitWindow(const double Dist,
			    std::vector<int>& window,
			 Geometry::Vec3D& Pt) const
  /*!
    Get the exit window
    \param Dist :: Dist from exit point
    \param window :: window vector of paired planes
    \param Pt :: Output point for tally
    \return Master Surface
   */
{
  ELog::RegMethod RegA("zoomInsertBlock","exitWindow");

  window.clear();
  window.push_back(SMap.realSurf(surfIndex+3));
  window.push_back(SMap.realSurf(surfIndex+4));
  window.push_back(SMap.realSurf(surfIndex+5));
  window.push_back(SMap.realSurf(surfIndex+6));

  Pt=Origin+Y*Dist;
  return SMap.realSurf(surfIndex+1);
}


} // NAMESPACE shutterSystem
