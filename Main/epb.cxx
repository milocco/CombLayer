/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   Main/epb.cxx
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
#include "MersenneTwister.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "InputControl.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Tensor.h"
#include "Vec3D.h"
#include "inputParam.h"
#include "Triple.h"
#include "NRange.h"
#include "NList.h"
#include "Tally.h"
#include "TallyCreate.h"
#include "Transform.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "surfIndex.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "ModeCard.h"
#include "PhysCard.h"
#include "LSwitchCard.h"
#include "PhysImp.h"
#include "KGroup.h"
#include "Source.h"
#include "KCode.h"
#include "PhysicsCards.h"
#include "BasicWWE.h"
#include "MainProcess.h"
#include "SimProcess.h"
#include "SimInput.h"
#include "SurInter.h"
#include "Simulation.h"
#include "SimPHITS.h"
#include "PointWeights.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinearComp.h"
#include "mainJobs.h"
#include "Volumes.h"
#include "DefPhysics.h"
#include "variableSetup.h"
#include "ImportControl.h"
#include "SourceCreate.h"
#include "SourceSelector.h"
#include "TallySelector.h"
#include "World.h"
#include "makeEPB.h"

MTRand RNG(12345UL);

///\cond STATIC
namespace ELog 
{
  ELog::OutputLog<EReport> EM;
  ELog::OutputLog<FileReport> FM("Spectrum.log");
  ELog::OutputLog<FileReport> RN("Renumber.txt");   ///< Renumber
  ELog::OutputLog<StreamReport> CellM;
}
///\endcond STATIC

int 
main(int argc,char* argv[])
{
  int exitFlag(0);                // Value on exit
  ELog::RegMethod RControl("","main");
  mainSystem::activateLogging(RControl);

  std::string Oname;
  std::vector<std::string> Names;  
  std::map<std::string,std::string> Values;  
  std::map<std::string,std::string> AddValues;  
  std::map<std::string,double> IterVal;           // Variable to iterate 

  // PROCESS INPUT:
  InputControl::mainVector(argc,argv,Names);
  mainSystem::inputParam IParam;
  createEPBInputs(IParam);

  const int iteractive(IterVal.empty() ? 0 : 1);   
  Simulation* SimPtr=createSimulation(IParam,Names,Oname);
  if (!SimPtr) return -1;

  // The big variable setting
  setVariable::EPBVariables(SimPtr->getDataBase());
  InputModifications(SimPtr,IParam,Names);
  
  // Definitions section 
  int MCIndex(0);
  const int multi=IParam.getValue<int>("multi");
  try
    {
      while(MCIndex<multi)
	{
	  if (MCIndex)
	    {
	      ELog::EM.setActive(4);    // write error only
	      ELog::FM.setActive(4);    
	      ELog::RN.setActive(0);    
	      // if (iteractive)
	      // 	mainSystem::incRunTimeVariable
	      // 	  (SimPtr->getDataBase(),IterVal);
	    }

	  SimPtr->resetAll();

	  epbSystem::makeEPB EPBObj;
	  World::createOuterObjects(*SimPtr);
	  EPBObj.build(SimPtr,IParam);

	  SDef::sourceSelection(*SimPtr,IParam);

	  SimPtr->removeComplements();
	  SimPtr->removeDeadSurfaces(0);         

	  ModelSupport::setDefaultPhysics(*SimPtr,IParam);
	  const int renumCellWork=tallySelection(*SimPtr,IParam);
	  SimPtr->masterRotation();
	  if (createVTK(IParam,SimPtr,Oname))
	    {
	      delete SimPtr;
	      ModelSupport::objectRegister::Instance().reset();
	      return 0;
	    }
	  if (IParam.flag("endf"))
	    SimPtr->setENDF7();
	  createMeshTally(IParam,SimPtr);

	  SimProcess::importanceSim(*SimPtr,IParam);
	  SimProcess::inputPatternSim(*SimPtr,IParam); // energy cut etc

	  if (renumCellWork)
	    tallyRenumberWork(*SimPtr,IParam);
	  tallyModification(*SimPtr,IParam);

	  if (IParam.flag("cinder"))
	    SimPtr->setForCinder();

	  // Ensure we done loop
	  do
	    {
	      SimProcess::writeIndexSim(*SimPtr,Oname,MCIndex);
	      MCIndex++;
	    }
	  while(!iteractive && MCIndex<multi);
	}
      if (IParam.flag("cinder"))
	SimPtr->writeCinder();
      ModelSupport::calcVolumes(SimPtr,IParam);
      ModelSupport::objectRegister::Instance().write("ObjectRegister.txt");
    }
  catch (ColErr::ExitAbort& EA)
    {
      if (!EA.pathFlag())
	ELog::EM<<"Exiting from "<<EA.what()<<ELog::endCrit;
      exitFlag=-2;
    }
  catch (ColErr::ExBase& A)
    {
      ELog::EM<<"EXCEPTION FAILURE :: "
	      <<A.what()<<ELog::endCrit;
      exitFlag= -1;
    }
  delete SimPtr;
  ModelSupport::objectRegister::Instance().reset();
  ModelSupport::surfIndex::Instance().reset();
  return exitFlag;
}
