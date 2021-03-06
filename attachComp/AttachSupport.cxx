/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   attachComp/AttachSupport.cxx
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
#include <deque>
#include <string>
#include <algorithm>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "Triple.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Surface.h"
#include "Rules.h"
#include "HeadRule.h"
#include "Object.h"
#include "surfRegister.h"
#include "objectRegister.h"

#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "NRange.h"
#include "NList.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "objectRegister.h"
#include "Qhull.h"
#include "KGroup.h"
#include "Source.h"
#include "Simulation.h"
#include "SurInter.h"
#include "AttachSupport.h"

#include "Debug.h"

namespace attachSystem
{

void createAddition(const int,Rule*,Rule*&);

void
addUnion(const MonteCarlo::Object& Obj,Rule*& outRule)
  /*!
    Given an Rule: intersect into it another object
    \param Obj :: Object add 
    \param outRule :: Rule to modify
  */
{
  ELog::RegMethod RegA("AttachSupport[F]","addUnion<Obj,Rule>");

  createAddition(-1,Obj.topRule()->clone(),outRule);
  return;
}

void
addUnion(const int SN,const Geometry::Surface* SPtr,
	 Rule*& outRule)
  /*!
    Given an Rule: intersect into it another object
    \param Obj :: Object add 
    \param outRule :: Rule to modify
  */
{
  ELog::RegMethod RegA("AttachSupport[F]","addUnion<int,Rule>");
  createAddition(-1,new SurfPoint(SPtr,SN),outRule);  
  return;
}

void
addIntersection(const MonteCarlo::Object& Obj,Rule*& outRule)
  /*!
    Given an Rule: intersect into it another object
    \param Obj :: Object add 
    \param outRule :: Rule to modify
  */
{
  ELog::RegMethod RegA("AttachSupport[F]","addIntersection<Obj,Rule>");

  createAddition(1,Obj.topRule()->clone(),outRule);
  return;
}

void
addIntersection(const int SN,const Geometry::Surface* SPtr,
		Rule*& outRule)
  /*!
    Given an Rule: intersect into it another object
    \param SN :: Surface number [signed]
    \param SPtr :: Surface Pointer
    \param outRule :: Rule to modify
  */
{
  ELog::RegMethod RegA("AttachSupport[F]","addIntersection<int,Rule>");
  createAddition(1,new SurfPoint(SPtr,SN),outRule);  
  return;
}

void
createAddition(const int InterFlag,Rule* NRptr,
	       Rule*& outRule)
  /*!
    Function to actually do the addition of a rule to 
    avoid code repeat.
    \param InterFlag :: Intersection / Union [ 1 : -1 ] : 0 for new 
    \param NRPtr :: New Rule pointer to add
    \param outRule :: Rule system to modify
   */
{
  ELog::RegMethod RegA("AttachSupport","createAddition");

  // Null object : Set and return
  if (!outRule || !InterFlag)
    {
      delete outRule;
      outRule=NRptr;
      return;
    }
  // This is an intersection and we want to add our rule at the base
  // Find first item that is not an intersection
  Rule* RPtr(outRule);
  
  std::deque<Rule*> curLevel;
  curLevel.push_back(outRule);
  while(!curLevel.empty())
    {
      RPtr=curLevel.front();
      curLevel.pop_front();
      if (RPtr->type() != InterFlag)  // Success not an intersection
	{
	  Rule* parent=RPtr->getParent();     // grandparent
	  Rule* sA=RPtr;
	  Rule* Item(0);
	  if (InterFlag==1)
	    Item=new Intersection(parent,sA,NRptr);
	  else
	    Item=new Union(parent,sA,NRptr);
	  // Find place ot insert it
	  if (!parent)
	    outRule=Item;
	  else
	    parent->setLeaf(Item,parent->findLeaf(RPtr));

	  return;
	}
      Rule* LeafPtr;
      if ( (LeafPtr=RPtr->leaf(0)) )
	curLevel.push_back(LeafPtr);
      if ( (LeafPtr=RPtr->leaf(1)) )
	curLevel.push_back(LeafPtr);
    }
  ELog::EM<<"Failed on rule structure  "<<ELog::endErr;
  return;  
}

void
addToInsertControl(Simulation& System,
		   const attachSystem::FixedComp& BaseFC,
		   const attachSystem::FixedComp& FC,
		   attachSystem::ContainedComp& CC)
/*!
  Adds this object to the containedComp to be inserted.
  FC is the fixed object that is to be inserted -- linkpoints
  must be set. It is tested against all the ojbect with
  this object .
  \param System :: Simulation to use
  \param FC :: FixedComp with the points
  \param CC :: ContainedComp object to add to this
*/
{
  ELog::RegMethod RegA("AttachSupport","addToInsertControl");
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  const int cellN=OR.getCell(BaseFC.getKeyName());
  const int cellR=OR.getRange(BaseFC.getKeyName());
  addToInsertControl(System,cellN,cellN+cellR,FC,CC);

  return;
}

void
addToInsertControl(Simulation& System,
		   const std::string& OName,
		   const attachSystem::FixedComp& FC,
		   attachSystem::ContainedComp& CC)
/*!
  Adds this object to the containedComp to be inserted.
  FC is the fixed object that is to be inserted -- linkpoints
  must be set. It is tested against all the ojbect with
  this object .
  \param System :: Simulation to use
  \param CellA :: First cell number
  \param CellB :: Last cell number
  \param FC :: FixedComp with the points
  \param CC :: ContainedComp object to add to this
*/
{
  ELog::RegMethod RegA("AttachSupport","addToInsertControl");
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  const int cellN=OR.getCell(OName);
  addToInsertControl(System,cellN,cellN+OR.getRange(OName),FC,CC);

  return;
}

void
addToInsertControl(Simulation& System,
		   const int cellA,const int cellB,
		   const attachSystem::FixedComp& FC,
		   attachSystem::ContainedComp& CC)
/*!
    Adds this object to the containedComp to be inserted.
    FC is the fixed object that is to be inserted -- linkpoints
    must be set. It is tested against all the ojbect with
    this object .
    \param System :: Simulation to use
    \param CellA :: First cell number [to test]
    \param CellB :: Last cell number  [to test]
    \param FC :: FixedComp with the points
    \param CC :: ContainedComp object to add to this
  */
{
  ELog::RegMethod RegA("AttachSupport","addToInsertControl");

  const size_t NPoint=FC.NConnect();
  for(int i=cellA+1;i<=cellB;i++)
    {
      MonteCarlo::Qhull* CRPtr=System.findQhull(i);
      if (i==cellA+1 && !CRPtr)
	throw ColErr::InContainerError<int>(i,"Object not build");
      else if (!CRPtr)
	break;

      CRPtr->populate();
      for(size_t j=0;j<NPoint;j++)
	{
	  const Geometry::Vec3D& Pt=FC.getLinkPt(j);
	  if (CRPtr->isValid(Pt))
	    {
	      CC.addInsertCell(i);
	      break;
	    }
	}
    }
  CC.insertObjects(System);
  return;
}


// SURFACE INTERSECT
void
addToInsertSurfCtrl(Simulation& System,
		    const attachSystem::FixedComp& BaseFC,
		    attachSystem::ContainedComp& CC)
  /*!
    Adds this object to the containedComp to be inserted.
    FC is the fixed object that is to be inserted -- linkpoints
    must be set. It is tested against all the ojbect with
    this object .
    \param System :: Simulation to use
    \param BaseFC :: FixedComp for name
    \param CC :: ContainedComp object to add to this
  */
{
  ELog::RegMethod RegA("AttachSupport","addToInsertSurfCtrl(FC,CC)");
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  const int cellN=OR.getCell(BaseFC.getKeyName());
  const int cellR=OR.getRange(BaseFC.getKeyName());
  addToInsertSurfCtrl(System,cellN,cellN+cellR,CC);

  return;
}

void
addToInsertSurfCtrl(Simulation& System,
		    const int cellA,const int cellB,
		    attachSystem::ContainedComp& CC)
 /*!
   Adds this object to the containedComp to be inserted.
   FC is the fixed object that is to be inserted -- linkpoints
   must be set. It is tested against all the ojbect with
   this object .
   \param System :: Simulation to use
   \param CellA :: First cell number [to test]
   \param CellB :: Last cell number  [to test]
   \param CC :: ContainedComp object to add to this
  */
{
  ELog::RegMethod RegA("AttachSupport","addToInsertSurfCtrl(int,int,CC)");

  const std::vector<Geometry::Surface*> SVec=CC.getSurfaces();

  std::vector<Geometry::Vec3D> Out;		
  std::vector<Geometry::Vec3D>::const_iterator vc;

  for(int i=cellA+1;i<=cellB;i++)
    {
      MonteCarlo::Qhull* CRPtr=System.findQhull(i);
      if (i==cellA+1 && !CRPtr)
	throw ColErr::InContainerError<int>(i,"Object not build");
      else if (!CRPtr)
	break;
      
      CRPtr->populate();
      CRPtr->createSurfaceList();
      const std::vector<const Geometry::Surface*>&
	CellSVec=CRPtr->getSurList();
      // LOOP OVER ALL SURFACE SETS:

      bool insertFlag(1);
      for(size_t iA=0;insertFlag && iA<SVec.size();iA++)
	for(size_t iB=0;insertFlag && iB<CellSVec.size();iB++)
	  for(size_t iC=iB+1;insertFlag && iC<CellSVec.size();iC++)
	    {
	      Out=SurInter::processPoint(SVec[iA],CellSVec[iB],CellSVec[iC]);
	      for(vc=Out.begin();vc!=Out.end();vc++)
		{
		  std::set<int> boundarySet;
		  boundarySet.insert(CellSVec[iB]->getName());
		  boundarySet.insert(CellSVec[iC]->getName());
		  if (CRPtr->isValid(*vc,boundarySet) &&
		      !CC.isOuterValid(*vc,SVec[iA]->getName()))
		    {
		      CC.addInsertCell(i);
		      insertFlag=0;
		      break;
		    }
		}
	    }
      if (insertFlag)  // No match found check link points:
	{
	  for(size_t iA=0;insertFlag && iA<SVec.size();iA++)
	    for(size_t iB=iA+1;insertFlag && iB<SVec.size();iB++)
	      for(size_t iC=iB+1;insertFlag && iC<SVec.size();iC++)
		{
		  Out=SurInter::processPoint(SVec[iA],SVec[iB],SVec[iC]);
		  for(vc=Out.begin();vc!=Out.end();vc++)
		    {
		      if (CRPtr->isValid(*vc))
			{
			  CC.addInsertCell(i);
			  insertFlag=0;
			  break;
			}
		    }
		}
	}
      if (insertFlag)  // No match found Now inter chekc
	{
	  for(size_t iA=0;insertFlag && iA<SVec.size();iA++)
	    for(size_t iB=iA+1;insertFlag && iB<SVec.size();iB++)
	      for(size_t iC=0;insertFlag && iC<CellSVec.size();iC++)
		{
		  Out=SurInter::processPoint(SVec[iA],SVec[iB],CellSVec[iC]);
		  for(vc=Out.begin();vc!=Out.end();vc++)
		    {
		      std::set<int> boundarySet;
		      boundarySet.insert(SVec[iA]->getName());
		      boundarySet.insert(SVec[iB]->getName());
		      if (CRPtr->isValid(*vc,CellSVec[iC]->getName()) &&
			  !CC.isOuterValid(*vc,boundarySet))
			{
			  CC.addInsertCell(i);
			  insertFlag=0;
			  break;
			}
		    }
		}
	}
    }
      
  CC.insertObjects(System);
  return;
}

void
addToInsertForced(Simulation& System,
		  const int cellA,const int cellB,
		  attachSystem::ContainedComp& CC)
 /*!
   Force CC into the BaseFC objects
  \param System :: Simulation to use
  \param BaseFC :: FixedComp object to have CC inserted into it.
  \param CC :: ContainedComp object to add to the BaseFC
 */
{
  ELog::RegMethod RegA("AttachSupport","addToInsertForce(int,int,CC)");

  for(int i=cellA+1;i<=cellB;i++)
    {
      MonteCarlo::Qhull* CRPtr=System.findQhull(i);
      if (i==cellA+1 && !CRPtr)
	throw ColErr::InContainerError<int>(i,"Object not built");
      else if (!CRPtr)
	break;
      CRPtr->populate();
      CRPtr->createSurfaceList();
      CC.addInsertCell(i);
    }

  CC.insertObjects(System);

  return;
}  


void
addToInsertForced(Simulation& System,
		  const attachSystem::FixedComp& BaseFC,
		  attachSystem::ContainedComp& CC)
 /*!
   Force CC into the BaseFC objects
  \param System :: Simulation to use
  \param CellA :: First cell number [to test]
   \param CellB :: Last cell number  [to test]
   
  \param CC :: ContainedComp object to add to the BaseFC
 */
{
  ELog::RegMethod RegA("AttachSupport","addToInsertForced(FC,CC)");
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  const int cellN=OR.getCell(BaseFC.getKeyName());
  const int cellR=OR.getRange(BaseFC.getKeyName());
  addToInsertForced(System,cellN,cellN+cellR,CC);
  return;
}  



}  // NAMESPACE attachSystem
