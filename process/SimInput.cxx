/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   process/SimInput.cxx
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
#include <map> 
#include <list> 
#include <set>
#include <string>
#include <algorithm>
#include <iterator>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

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
#include "inputParam.h"
#include "Triple.h"
#include "NList.h"
#include "NRange.h"
#include "Rules.h"
#include "Code.h"
#include "FItem.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "ModeCard.h"
#include "PhysCard.h"
#include "PhysImp.h"
#include "LSwitchCard.h"
#include "KGroup.h"
#include "Source.h"
#include "KCode.h"
#include "PhysicsCards.h"
#include "Simulation.h"
#include "ImportControl.h"
#include "SimValid.h"
#include "MainProcess.h"
#include "BasicWWE.h"
#include "SimInput.h"


namespace SimProcess
{

void
importanceSim(Simulation& System,const mainSystem::inputParam& IParam)
  /*!
    Apply importances/renumber and weights
     \param System :: Simuation object 
     \param IParam :: Inputer
   */
{
  ELog::RegMethod RegA("SimInput","importanceSim");

  System.populateCells();
  System.createObjSurfMap();

  WeightSystem::simulationImp(System,IParam);
  mainSystem::renumberCells(System,IParam);
  WeightSystem::simulationWeights(System,IParam);
  return;
}

void
inputPatternSim(Simulation& System,const mainSystem::inputParam& IParam)
  /*!
    Check the validity of the simulation
    \param System :: Simuation object 
    \param IParam :: Inpute parameters
  */
{
  ELog::RegMethod RegA("SimInput","inputPatterSim");
  if (IParam.flag("validCheck"))
    {
      ELog::EM<<"TRACK "<<ELog::endDebug;
      ModelSupport::SimValid SValidCheck;
      SValidCheck.run(System,IParam.getValue<size_t>("validCheck"));
    }
  
  // RENUMBER:
  if (IParam.flag("cinder"))
    System.setForCinder();
  
  // Cut energy tallies:
  if (IParam.flag("ECut"))
    System.setEnergy(IParam.getValue<double>("ECut"));
  
  return;
}

}  // NAMESPACE SimProcess
