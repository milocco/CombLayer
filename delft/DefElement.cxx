/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   delft/DefElement.cxx
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
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

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

#include "ReactorGrid.h"
#include "RElement.h"
#include "DefElement.h"

namespace delftSystem
{


DefElement::DefElement(const size_t XI,const size_t YI,
		       const std::string& Key) :
  RElement(XI,YI,Key)
  /*!
    Constructor BUT ALL variable are left unpopulated.
  */
{}


void
DefElement::populate(const Simulation& System)
  /*!
    Populate all the variables
    Requires that unset values are copied from previous block
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("DefElement","populate");

  return;
}

void
DefElement::createUnitVector(const FixedComp& FC,
			      const Geometry::Vec3D& OG)
  /*!
    Create the unit vectors
    - Y Down the beamline
    \param FC :: Reactor Grid Unit
    \param OG :: Orgin
  */
{
  ELog::RegMethod RegA("DefElement","createUnitVector");

  attachSystem::FixedComp::createUnitVector(FC);
  Origin=OG;
  return;
}

void
DefElement::createSurfaces()
  /*!
    Creates/duplicates the surfaces for this block
    \param startSurf :: First surface to use
  */
{  
  ELog::RegMethod RegA("ZoomInsertBlock","createSurface");


  return;
}

void
DefElement::createObjects(Simulation& System)
  /*!
    Create the objects
    \param System :: Simulation
  */
{
  ELog::RegMethod RegA("DefElement","createObjects");

  
  return;
}


void
DefElement::createLinks()
  /*!
    Creates a full attachment set
    0 - 1 standard points
    2 - 3 beamaxis points
  */
{

  return;
}

void
DefElement::createAll(Simulation& System,const FixedComp& FC,
		       const Geometry::Vec3D& OG)
  /*!
    Global creation of the hutch
    \param System :: Simulation to add vessel to
    \param FC :: Fixed Unit
    \param OG :: Orgin
  */
{
  ELog::RegMethod RegA("DefElement","createAll(DefElement)");
  populate(System);

  createUnitVector(FC,OG);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);       

  return;
}


} // NAMESPACE delftSystem
