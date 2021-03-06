/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuild/essVariables.cxx
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
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "support.h"
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "variableSetup.h"

namespace setVariable
{
  void EssBeamLinesVariables(FuncDataBase&);
  void EssConicModerator(FuncDataBase&);
  void ESSWheel(FuncDataBase&);
  void ESSLayerMod(FuncDataBase&);

void ESSLayerMod(FuncDataBase& Control)
  /*
    Set the variables for the lower moderators
    \param Control :: DataBase to put variables
   */
{
  Control.addVariable("LowModXStep",0.0);  
  Control.addVariable("LowModYStep",0.0);  
  Control.addVariable("LowModZStep",-18.0);
  Control.addVariable("LowModXYangle",125.15); 
  Control.addVariable("LowModZangle",0.0);

  Control.addVariable("LowModRadius",8.0);
  Control.addVariable("LowModHeight",8.0);
  Control.addVariable("LowModMat",1001);
  Control.addVariable("LowModTemp",300.0);

  Control.addVariable("LowModNLayers",6);
  // al layer
  Control.addVariable("LowModHGap1",0.3);
  Control.addVariable("LowModRadGap1",0.3);
  Control.addVariable("LowModMaterial1",13061);  // Al materk
  Control.addVariable("LowModTemp1",20.0);  
  // Vac gap
  Control.addVariable("LowModHGap2",0.5);
  Control.addVariable("LowModRadGap2",0.5);
  Control.addVariable("LowModMaterial2",0); 
  Control.addVariable("LowModTemp3",300.0);
  // Next Al layer
  Control.addVariable("LowModHGap3",0.2);
  Control.addVariable("LowModRadGap3",0.2);
  Control.addVariable("LowModMaterial3",13060); 
  Control.addVariable("LowModTemp3",300.0);  
  // He Layer
  Control.addVariable("LowModHGap4",0.2);
  Control.addVariable("LowModRadGap4",0.2);
  Control.addVariable("LowModMaterial4",0); 
  // Outer Layer
  Control.addVariable("LowModHGap5",0.2);
  Control.addVariable("LowModRadGap5",0.2);
  Control.addVariable("LowModMaterial5",5); 
  Control.addVariable("LowModTemp5",300.0); 
  // Clearance
  Control.addVariable("LowModHGap6",0.2);
  Control.addVariable("LowModRadGap6",0.2);
  Control.addVariable("LowModMaterial6",0); 
  return;
}



void
EssWheel(FuncDataBase& Control)
  /*!
    Variables that are used for the segmented wheel
    \param Control :: Segment variables
   */
{
  // WHEEL SHAFT

  Control.addVariable("SegWheelShaftNLayers",6);
  Control.addVariable("SegWheelShaftTopHeight",400.0);
  Control.addVariable("SegWheelShaftBaseHeight",-13.5);
  // Control.addVariable("SegWheelShaftHeight",435.0);

  // Control.addVariable("SegWheelShaftRadius",25.0);
  // see WheelInnerRadius below (28cm): it's the external joint radius
  Control.addVariable("SegWheelShaftJointThick",3.0);
  Control.addVariable("SegWheelShaftCoolThick",4.0);
  Control.addVariable("SegWheelShaftCladThick",3.0);
  Control.addVariable("SegWheelShaftVoidThick",1.0);
  Control.addVariable("SegWheelShaftSupportRadius",40.0);
  Control.addVariable("SegWheelShaftSupportThick",2.5);
  Control.addVariable("SegWheelShaftBaseThick",3.0);
  Control.addVariable("SegWheelShaftBaseFootThick",13.5);

  Control.addVariable("SegWheelCladShaftMat",26317);
  Control.addVariable("SegWheelCoolingShaftMatInt",2002);
  Control.addVariable("SegWheelCoolingShaftMatExt",2000);
 
  // TARGET

  Control.addVariable("SegWheelXStep",0.0);  
  Control.addVariable("SegWheelYStep",115.0);  
  Control.addVariable("SegWheelZStep",0.0);
  Control.addVariable("SegWheelXYangle",0.0); 
  Control.addVariable("SegWheelZangle",0.0);
  //
  Control.addVariable("SegWheelTargetHeight",8.0);
  Control.addVariable("SegWheelTargetSectorOffsetX",0.0);  
  Control.addVariable("SegWheelTargetSectorOffsetY",142.6476749356);
  Control.addVariable("SegWheelTargetSectorOffsetZ",0.0);
  Control.addVariable("SegWheelTargetSectorAngleXY",1.1234180949);
  Control.addVariable("SegWheelTargetSectorAngleZ",0.0);
  Control.addVariable("SegWheelTargetSectorApertureXY",2.2468361899);
  Control.addVariable("SegWheelTargetSectorApertureZ",0.0);
  Control.addVariable("SegWheelTargetSectorNumber",33);
  //
  Control.addVariable("SegWheelCoolantThickOut",0.15);
  Control.addVariable("SegWheelCoolantThickIn",0.65);
  //
  Control.addVariable("SegWheelCaseThickZ",0.5);
  Control.addVariable("SegWheelCaseThickX",0.1);
  //
  Control.addVariable("SegWheelVoidThick",2.15);  //distance target-ion tube height

  Control.addVariable("SegWheelInnerRadius",28.0);
  Control.addVariable("SegWheelCoolantRadiusOut",124.8);
  Control.addVariable("SegWheelCoolantRadiusIn",114.5);
  Control.addVariable("SegWheelCaseRadius",125.0);
  Control.addVariable("SegWheelVoidRadius",126.0);

  Control.addVariable("SegWheelWMat",74001);
  Control.addVariable("SegWheelSteelMat",26317);
  Control.addVariable("SegWheelHeMat",2002);

  Control.addVariable("SegWheelInnerMat",3);

  Control.addVariable("SegWheelNLayers",25);

  Control.addVariable("SegWheelRadius1",84.27);
  Control.addVariable("SegWheelMatTYPE1",1);

  Control.addVariable("SegWheelRadius2",85.65);
  Control.addVariable("SegWheelMatTYPE2",2);

  Control.addVariable("SegWheelRadius3",97.02);
  Control.addVariable("SegWheelMatTYPE3",3);

  Control.addVariable("SegWheelRadius4",97.52);
  Control.addVariable("SegWheelMatTYPE4",2);

  Control.addVariable("SegWheelRadius5",102.37); 
  Control.addVariable("SegWheelMatTYPE5",3);

  Control.addVariable("SegWheelRadius6",102.77);
  Control.addVariable("SegWheelMatTYPE6",2);

  Control.addVariable("SegWheelRadius7",106.97);
  Control.addVariable("SegWheelMatTYPE7",3);

  Control.addVariable("SegWheelRadius8",107.29);
  Control.addVariable("SegWheelMatTYPE8",2);

  Control.addVariable("SegWheelRadius9",110.49);
  Control.addVariable("SegWheelMatTYPE9",3);

  Control.addVariable("SegWheelRadius10",110.78);
  Control.addVariable("SegWheelMatTYPE10",2); 

  Control.addVariable("SegWheelRadius11",112.78);
  Control.addVariable("SegWheelMatTYPE11",3); 

  Control.addVariable("SegWheelRadius12",113.05);
  Control.addVariable("SegWheelMatTYPE12",2); 

  Control.addVariable("SegWheelRadius13",114.65);
  Control.addVariable("SegWheelMatTYPE13",3); 

  Control.addVariable("SegWheelRadius14",114.9);
  Control.addVariable("SegWheelMatTYPE14",2); 

  // Control.addVariable("SegWheelRadius15",114.9);
  // Control.addVariable("SegWheelMatTYPE15",1); 

  Control.addVariable("SegWheelRadius15",116.5);
  Control.addVariable("SegWheelMatTYPE15",3); 

  Control.addVariable("SegWheelRadius16",116.75);
  Control.addVariable("SegWheelMatTYPE16",2); 

  Control.addVariable("SegWheelRadius17",118.35);
  Control.addVariable("SegWheelMatTYPE17",3); 

  Control.addVariable("SegWheelRadius18",118.6);
  Control.addVariable("SegWheelMatTYPE18",2); 

  Control.addVariable("SegWheelRadius19",120.0);
  Control.addVariable("SegWheelMatTYPE19",3); 

  Control.addVariable("SegWheelRadius20",120.25);
  Control.addVariable("SegWheelMatTYPE20",2); 

  Control.addVariable("SegWheelRadius21",121.65);
  Control.addVariable("SegWheelMatTYPE21",3); 
 // SHUTTER BAY
  Control.addVariable("ShutterBayXStep",0.0);  
  Control.addVariable("ShutterBayYStep",0.0);  
  Control.addVariable("ShutterBayZStep",0.0);
  Control.addVariable("ShutterBayXYangle",0.0); 
  Control.addVariable("ShutterBayZangle",0.0);
  Control.addVariable("ShutterBayRadius",600.0);
  Control.addVariable("ShutterBayHeight",400.0);
  Control.addVariable("ShutterBayDepth",400.0);
  Control.addVariable("ShutterBayMat",26000);

  Control.addVariable("SegWheelRadius22",121.9);
  Control.addVariable("SegWheelMatTYPE22",2); 

  Control.addVariable("SegWheelRadius23",123.1);
  Control.addVariable("SegWheelMatTYPE23",3); 

  Control.addVariable("SegWheelRadius24",123.35);
  Control.addVariable("SegWheelMatTYPE24",2); 

  Control.addVariable("SegWheelRadius25",124.55);
  Control.addVariable("SegWheelMatTYPE25",3); 
  return;
}


void
EssVariables(FuncDataBase& Control)
  /*!
    Function to set the control variables and constants
    -- This version is for ESS ()
    \param Control :: Function data base to add constants too
  */
{
// -----------
// GLOBAL stuff
// -----------

  Control.addVariable("zero",0.0);     // Zero
  Control.addVariable("one",1.0);      // one

 ///////////////////////////////////////


  // LOW MODERATOR
  Control.addVariable("LowModXStep",0.0);  
  Control.addVariable("LowModYStep",0.0);  
  Control.addVariable("LowModZStep",-18.0);
  Control.addVariable("LowModXYangle",125.15); 
  Control.addVariable("LowModZangle",0.0);
  Control.addVariable("LowModRadius",8.0);
  Control.addVariable("LowModHeight",13.0);
  Control.addVariable("LowModMat",1001);
  Control.addVariable("LowModTemp",300.0);

  Control.addVariable("LowModNLayers",3);
  // al layer
  Control.addVariable("LowModHGap1",0.3);
  Control.addVariable("LowModRadGap1",0.3);
  Control.addVariable("LowModMaterial1",13061);  // Al materk
  Control.addVariable("LowModTemp1",300.0);  
  // Vac gap
  Control.addVariable("LowModHGap2",0.5);
  Control.addVariable("LowModRadGap2",0.5);
  Control.addVariable("LowModMaterial2",0); 
 Control.addVariable("LowModTemp2",300.0);  
 // Outer Layer
  Control.addVariable("LowModHGap3",0.2);
  Control.addVariable("LowModRadGap3",0.2);
  Control.addVariable("LowModMaterial3",13060); 
  Control.addVariable("LowModTemp3",300.0);  

  // LOW PREMODERATOR 
  Control.addVariable("LowPreNLayers",3);  
  Control.addVariable("LowPreHeight1",2.0);  
  Control.addVariable("LowPreDepth1",1.0);  
  Control.addVariable("LowPreThick1",1.0);  
  Control.addVariable("LowPreMaterial1",1015);  
  Control.addVariable("LowPreHeight2",0.2);  
  Control.addVariable("LowPreDepth2",0.2);  
  Control.addVariable("LowPreThick2",0.2);  
  Control.addVariable("LowPreMaterial2",13060); 
  Control.addVariable("LowPreHeight3",2.15);  
  Control.addVariable("LowPreDepth3",0.0);  
  Control.addVariable("LowPreThick3",0.0);  
  Control.addVariable("LowPreMaterial3",2000);  
  // Control.addVariable("LowPreHeight4",0.0);  
  // Control.addVariable("LowPreDepth4",2.15);  
  // Control.addVariable("LowPreThick4",0.0);  
  // Control.addVariable("LowPreMaterial4",2000);  
 
  Control.addVariable("LowPreNView",2);  
  Control.addVariable("LowPreViewHeight1",12.0);  
  Control.addVariable("LowPreViewWidth1",13.8);  
  Control.addVariable("LowPreViewAngle1",0.0);  
  Control.addVariable("LowPreViewOpenAngle1",30.0);  

  Control.addVariable("LowPreViewHeight2",12.0);  
  Control.addVariable("LowPreViewWidth2",13.85);  
  Control.addVariable("LowPreViewAngle2",175.2);  
  Control.addVariable("LowPreViewOpenAngle2",30.0);  

  // LOW THERMAL BLOCKS
  Control.addVariable("LowPreABlockSide",0);  
  Control.addVariable("LowPreABlockActive",1);  
  Control.addVariable("LowPreABlockWidth",4.0);  
  Control.addVariable("LowPreABlockHeight",18.0);  
  Control.addVariable("LowPreABlockLength",10.4);  
  Control.addVariable("LowPreABlockWallThick",0.2);  
  Control.addVariable("LowPreABlockGap",0.0);  
  Control.addVariable("LowPreABlockWallMat",13060);  
  Control.addVariable("LowPreABlockWaterMat",1015);  
  Control.addVariable("LowPreABlockPipeLength",65);  //Be reflector radius
  Control.addVariable("LowPreABlockPipeInnerRadius",1.3);
  Control.addVariable("LowPreABlockPipeN",3);

  Control.addVariable("LowPreBBlockSide",0);  
  Control.addVariable("LowPreBBlockActive",1);  
  Control.addVariable("LowPreBBlockWidth",4.0);  
  Control.addVariable("LowPreBBlockHeight",18.0);  
  Control.addVariable("LowPreBBlockLength",10.4);  
  Control.addVariable("LowPreBBlockWallThick",0.2);  
  Control.addVariable("LowPreBBlockGap",0.0);  
  Control.addVariable("LowPreBBlockWallMat",13060);  
  Control.addVariable("LowPreBBlockWaterMat",1015);  
  Control.addVariable("LowPreBBlockPipeLength",65);  //Be reflector radius
  Control.addVariable("LowPreBBlockPipeInnerRadius",1.3);
  Control.addVariable("LowPreBBlockPipeN",3);

 // LOW FLIGHTLINES
  Control.addVariable("LowAFlightXStep",-5.3);      // Step from centre
  Control.addVariable("LowAFlightZStep",0.0);      // Step from centre
  Control.addVariable("LowAFlightAngleXY1",26.5);  // Angle out
  Control.addVariable("LowAFlightAngleXY2",29.8);  // Angle out
  Control.addVariable("LowAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("LowAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("LowAFlightHeight",12.0);     // Full height
  Control.addVariable("LowAFlightWidth",25.3);     // Full width
  Control.addVariable("LowAFlightMaterial",2000);     // Helium

  Control.addVariable("LowBFlightXStep",6.55);     // Angle
  Control.addVariable("LowBFlightZStep",0.0);      // Step from centre
  Control.addVariable("LowBFlightAngleXY1",24.4);  // Angle out
  Control.addVariable("LowBFlightAngleXY2",32.0);  // Angle out
  Control.addVariable("LowBFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("LowBFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("LowBFlightHeight",12.0);     // Full height
  Control.addVariable("LowBFlightWidth",25.5);     // Full width
  Control.addVariable("LowBFlightMaterial",2000);     // Helium

 // LOW H2 PIPE
  Control.addVariable("LSupplyNSegIn",2);
  Control.addVariable("LSupplyPPt0",Geometry::Vec3D(0,0,0.0));
  Control.addVariable("LSupplyPPt1",Geometry::Vec3D(0,-19.25,0));
  // mod xyangle (125.5)+pipe xyangle (57.5)
  Control.addVariable("LSupplyPPt2",Geometry::Vec3D(-64.930,-19.25,3.005));
  Control.addVariable("LSupplyInRadius",1.5);
  Control.addVariable("LSupplyInAlRadius",1.7);
  Control.addVariable("LSupplyMidAlRadius",1.8);
  Control.addVariable("LSupplyVoidRadius",2.3);
  Control.addVariable("LSupplyOutAlRadius",2.5);
  Control.addVariable("LSupplyInMat",1001);
  Control.addVariable("LSupplyInAlMat",13061);
  Control.addVariable("LSupplyMidAlMat",13061);
  Control.addVariable("LSupplyVoidMat",0);
  Control.addVariable("LSupplyOutAlMat",13060);

  // low mod return pipe
  Control.addVariable("LReturnNSegIn",1);
  Control.addVariable("LReturnPPt0",Geometry::Vec3D(0.0,0.0,5.0));
  Control.addVariable("LReturnPPt1",Geometry::Vec3D(2.628,56.939,5.0));
  // Control.addVariable("LReturnPPt0",Geometry::Vec3D(-7.692,0.356,4.0));
  // Control.addVariable("LReturnPPt1",Geometry::Vec3D(-64.930,3.005,4.0));
  Control.addVariable("LReturnInRadius",1.5);
  Control.addVariable("LReturnInAlRadius",1.8);
  Control.addVariable("LReturnMidAlRadius",1.8);
  Control.addVariable("LReturnVoidRadius",2.3);
  Control.addVariable("LReturnOutAlRadius",2.5);
  Control.addVariable("LReturnInMat",1001);
  Control.addVariable("LReturnInAlMat",13061);
  Control.addVariable("LReturnMidAlMat",13061);
   Control.addVariable("LReturnVoidMat",0);
  Control.addVariable("LReturnOutAlMat",13060);


  // TOP MODERATOR
  Control.addVariable("TopModXStep",0.0);  
  Control.addVariable("TopModYStep",0.0);  
  Control.addVariable("TopModZStep",18.0);
  Control.addVariable("TopModXYangle",57.5); 
  Control.addVariable("TopModZangle",0.0);

  Control.addVariable("TopModRadius",8.0);
  Control.addVariable("TopModHeight",13.0);
  Control.addVariable("TopModMat",1001);
  Control.addVariable("TopModTemp",20.0);

  Control.addVariable("TopModNLayers",3);
  // al layer
  Control.addVariable("TopModHGap1",0.3);
  Control.addVariable("TopModRadGap1",0.3);
  Control.addVariable("TopModMaterial1",13061);  // Al materk
  Control.addVariable("TopModTemp1",0.0);  
  // Vac gap
  Control.addVariable("TopModHGap2",0.5);
  Control.addVariable("TopModRadGap2",0.5);
  Control.addVariable("TopModMaterial2",0); 
  // Next Al layer
  Control.addVariable("TopModHGap3",0.2);
  Control.addVariable("TopModRadGap3",0.2);
  Control.addVariable("TopModMaterial3",13060); 
  Control.addVariable("TopModTemp3",77.0);  

  // TOP PREMODERATOR
  Control.addVariable("TopPreNLayers",3);  
  Control.addVariable("TopPreHeight1",1.0);  
  Control.addVariable("TopPreDepth1",2.0);  
  Control.addVariable("TopPreThick1",1.0);  
  Control.addVariable("TopPreMaterial1",1015);  
  Control.addVariable("TopPreHeight2",0.2);  
  Control.addVariable("TopPreDepth2",0.2);  
  Control.addVariable("TopPreThick2",0.2);  
  Control.addVariable("TopPreMaterial2",13060);  
  Control.addVariable("TopPreHeight3",0.0);  
  Control.addVariable("TopPreDepth3",2.15);  
  Control.addVariable("TopPreThick3",0.0);  
  Control.addVariable("TopPreMaterial3",2000);  
 
  Control.addVariable("TopPreNView",2);  
  Control.addVariable("TopPreViewHeight1",12.0);  
  Control.addVariable("TopPreViewWidth1",14.2);  
  Control.addVariable("TopPreViewAngle1",2.75);  
  Control.addVariable("TopPreViewOpenAngle1",30.0);  

  Control.addVariable("TopPreViewHeight2",12.0);  
  Control.addVariable("TopPreViewWidth2",14.);  
  Control.addVariable("TopPreViewAngle2",177.25);  
  Control.addVariable("TopPreViewOpenAngle2",30.0);  

  // thermal blocks
  Control.addVariable("TopPreABlockSide",0);  
  Control.addVariable("TopPreABlockActive",1);  
  Control.addVariable("TopPreABlockWidth",4.0);  
  Control.addVariable("TopPreABlockHeight",18.0);  
  Control.addVariable("TopPreABlockLength",10.4);  
  Control.addVariable("TopPreABlockWallThick",0.2);  
  Control.addVariable("TopPreABlockGap",0.0);  
  Control.addVariable("TopPreABlockWallMat",13060);  
  Control.addVariable("TopPreABlockWaterMat",1015);
  Control.addVariable("TopPreABlockPipeLength",65);  //Be reflector radius
  Control.addVariable("TopPreABlockPipeInnerRadius",1.3);
  Control.addVariable("TopPreABlockPipeN",3);
  
  Control.addVariable("TopPreBBlockActive",1);  
  Control.addVariable("TopPreBBlockSide",0);  
  Control.addVariable("TopPreBBlockWidth",4.0);  
  Control.addVariable("TopPreBBlockHeight",18.0);  
  Control.addVariable("TopPreBBlockLength",10.4);  
  Control.addVariable("TopPreBBlockWallThick",0.2);  
  Control.addVariable("TopPreBBlockGap",0.0);  
  Control.addVariable("TopPreBBlockWallMat",13060);  
  Control.addVariable("TopPreBBlockWaterMat",1015);  
  Control.addVariable("TopPreBBlockPipeLength",65); //Be reflector radius
  Control.addVariable("TopPreBBlockPipeInnerRadius",1.3);
  Control.addVariable("TopPreBBlockPipeN",3);
  

  // TOP A FLIGHT
  Control.addVariable("TopAFlightXStep",-6);      // Step from centre
  Control.addVariable("TopAFlightZStep",0.0);      // Step from centre
  Control.addVariable("TopAFlightAngleXY1",29.22);  // Angle out
  Control.addVariable("TopAFlightAngleXY2",27.05);  // Angle out
  Control.addVariable("TopAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("TopAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("TopAFlightHeight",12.0);     // Full height
  Control.addVariable("TopAFlightWidth",25.8);     // Full width
  Control.addVariable("TopAFlightMaterial",2000);     // Helium
 

  // TOP B FLIGHT
  Control.addVariable("TopBFlightXStep",6.);     // Angle
  Control.addVariable("TopBFlightZStep",0.0);      // Step from centre
  Control.addVariable("TopBFlightAngleXY1",27.0);  // Angle out
  Control.addVariable("TopBFlightAngleXY2",29.25);  // Angle out
  Control.addVariable("TopBFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("TopBFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("TopBFlightHeight",12.0);     // Full height
  Control.addVariable("TopBFlightWidth",25.6);     // Full width
  Control.addVariable("TopBFlightMaterial",2000);     // Helium

  // TOP H2 SUPPLY
  Control.addVariable("TSupplyNSegIn",2);
  Control.addVariable("TSupplyPPt0",Geometry::Vec3D(0.0,0.0,0.0));
  Control.addVariable("TSupplyPPt1",Geometry::Vec3D(0.0,-19.25,0.0));
  // mod xyangle (125.5)+pipe xyangle (57.5)
  Control.addVariable("TSupplyPPt2",Geometry::Vec3D(-64.999,-19.25,0.0085));
  Control.addVariable("TSupplyInRadius",1.5);
  Control.addVariable("TSupplyInAlRadius",1.7);
  Control.addVariable("TSupplyMidAlRadius",1.8);
  Control.addVariable("TSupplyVoidRadius",2.3);
  Control.addVariable("TSupplyOutAlRadius",2.5);
  Control.addVariable("TSupplyInMat",1001);
  Control.addVariable("TSupplyInAlMat",13061);
  Control.addVariable("TSupplyMidAlMat",13061);
  Control.addVariable("TSupplyVoidMat",0);
  Control.addVariable("TSupplyOutAlMat",13060);

  // TOP H2 RETURN
  Control.addVariable("TReturnNSegIn",1);
  Control.addVariable("TReturnPPt0",Geometry::Vec3D(0.0,0.0,-5.0));
  Control.addVariable("TReturnPPt1",Geometry::Vec3D(0.0,57.0,-5.0));
  Control.addVariable("TReturnInRadius",1.5);
  Control.addVariable("TReturnInAlRadius",1.8);
  Control.addVariable("TReturnMidAlRadius",1.8);
  Control.addVariable("TReturnVoidRadius",2.3);
  Control.addVariable("TReturnOutAlRadius",2.5);
  Control.addVariable("TReturnInMat",1001);
  Control.addVariable("TReturnInAlMat",13061);
  Control.addVariable("TReturnMidAlMat",13061);
  Control.addVariable("TReturnVoidMat",0);
  Control.addVariable("TReturnOutAlMat",13060);


  // PROTON BEAM
  //  Control.addVariable("ProtonBeamViewRadius",4.0);  
  Control.addVariable("ProtonTubeXStep",0.0);  
  Control.addVariable("ProtonTubeYStep",0.0);  
  Control.addVariable("ProtonTubeZStep",0.0);
  Control.addVariable("ProtonTubeXYangle",0.0); 
  Control.addVariable("ProtonTubeZangle",0.0);

  Control.addVariable("ProtonTubeSectionN",4);

  Control.addVariable("ProtonTubeRadius1",11.5);
  Control.addVariable("ProtonTubeLength1",120.0); //from mod centre leftside
  Control.addVariable("ProtonTubeZcut1",5.35); //cut Z planes
  Control.addVariable("ProtonTubeWallThick1",0.0);
  Control.addVariable("ProtonTubeInnerMat1",2000);
  Control.addVariable("ProtonTubeWallMat1",26316);

  Control.addVariable("ProtonTubeRadius2",10.5);
  Control.addVariable("ProtonTubeLength2",200.0);
  Control.addVariable("ProtonTubeZcut2",0); 
  Control.addVariable("ProtonTubeWallThick2",1.0);
  Control.addVariable("ProtonTubeInnerMat2",2000);
  Control.addVariable("ProtonTubeWallMat2",26316);

  Control.addVariable("ProtonTubeRadius3",10.5);
  Control.addVariable("ProtonTubeLength3",127.5);
  Control.addVariable("ProtonTubeZcut3",0);
  Control.addVariable("ProtonTubeWallThick3",1.0);
  Control.addVariable("ProtonTubeInnerMat3",2000);
  Control.addVariable("ProtonTubeWallMat3",26316);

  Control.addVariable("ProtonTubeRadius4",10.5);
  // Control.addVariable("ProtonTubeLength4",147.5);
  Control.addVariable("ProtonTubeLength4",152.5);
  Control.addVariable("ProtonTubeZcut4",10.5);
  Control.addVariable("ProtonTubeWallThick4",1.0);
  Control.addVariable("ProtonTubeInnerMat4",0);
  Control.addVariable("ProtonTubeWallMat4",26316);


  Control.addVariable("ProtonTubeWindowOffsetX",0.0);
  Control.addVariable("ProtonTubeWindowOffsetY",-450.0);
  Control.addVariable("ProtonTubeWindowOffsetZ",0.0);
  
  Control.addVariable("ProtonTubeWindowBoxSide",70.0);
  Control.addVariable("ProtonTubeWindowBoxSectionsN",5); //other are symmetric

  Control.addVariable("ProtonTubeWindowBoxRadius1",11.5);
  Control.addVariable("ProtonTubeWindowBoxThick1",2.5);
  Control.addVariable("ProtonTubeWindowBoxMat1",2000);
  Control.addVariable("ProtonTubeWindowBoxRadius2",32.5);
  Control.addVariable("ProtonTubeWindowBoxThick2",10.0);
  Control.addVariable("ProtonTubeWindowBoxMat2",26000);
  Control.addVariable("ProtonTubeWindowBoxRadius3",20.0);
  Control.addVariable("ProtonTubeWindowBoxThick3",2.5);
  Control.addVariable("ProtonTubeWindowBoxMat3",26316);
  Control.addVariable("ProtonTubeWindowBoxRadius4",17.5);
  Control.addVariable("ProtonTubeWindowBoxThick4",17.5);
  Control.addVariable("ProtonTubeWindowBoxMat4",26316);
  Control.addVariable("ProtonTubeWindowBoxRadius5",25.0);
  Control.addVariable("ProtonTubeWindowBoxThick5",5.0);
  Control.addVariable("ProtonTubeWindowBoxMat5",13060);

  Control.addVariable("ProtonTubeWindowBoxSide5",50.0);
  Control.addVariable("ProtonTubeWindowBoxHeightA5",14.60);
  Control.addVariable("ProtonTubeWindowBoxHeightB5",8.60);
  Control.addVariable("ProtonTubeWindowBoxHeightC5",7.5);
  Control.addVariable("ProtonTubeWindowBoxHeightD5",7.0);

  Control.addVariable("ProtonTubeWindowBoxWidthA5",21.40);
  Control.addVariable("ProtonTubeWindowBoxWidthB5",20.30);
  Control.addVariable("ProtonTubeWindowBoxWidthC5",19.80);

  Control.addVariable("ProtonTubeWindowBoxThickA5",3.2);
  Control.addVariable("ProtonTubeWindowBoxThickB5",0.5);
  Control.addVariable("ProtonTubeWindowBoxThickC5",0.2);

  Control.addVariable("ProtonTubeWindowBoxTubeN",33);
  Control.addVariable("ProtonTubeWindowBoxTubeRadius",0.27);
  Control.addVariable("ProtonTubeWindowBoxTubeThick",0.03);

  // Control.addVariable("ProtonTubeWindowBoxTubeThick",0.03);
  Control.addVariable("ProtonTubeWindowBoxTubeHeMat",2001);
  Control.addVariable("ProtonTubeWindowBoxTubeAlMat",13060);
  Control.addVariable("ProtonTubeWindowBoxExtHeMat",2000);
 

  // Control.addVariable("ProtonTubeHeight",600.0);
  // Control.addVariable("ProtonTubeWallThick",0.0);
  // Control.addVariable("ProtonTubeRefMat",77);
  // Control.addVariable("ProtonTubeWallMat",77);
  // Control.addVariable("ProtonTubeZcut",-63.2455);

 

  // WHEEL SHAFT
  Control.addVariable("WheelShaftNLayers",3);
  Control.addVariable("WheelShaftHeight",435.0);
  Control.addVariable("WheelShaftRadius",32.0);
  Control.addVariable("WheelShaftCoolThick",1.0);
  Control.addVariable("WheelShaftCladThick",0.5);
  Control.addVariable("WheelShaftVoidThick",0.8);

  Control.addVariable("WheelCladShaftMat",26317);
  Control.addVariable("WheelMainShaftMat",2002);

  // WHEEL BODY
  Control.addVariable("WheelXStep",0.0);  
  Control.addVariable("WheelYStep",115.0);  
  Control.addVariable("WheelZStep",0.0);
  Control.addVariable("WheelXYangle",0.0); 
  Control.addVariable("WheelZangle",0.0);
  Control.addVariable("WheelTargetHeight",8.0);
  Control.addVariable("WheelCoolantThick",0.5);
  Control.addVariable("WheelCaseThick",0.5);
  Control.addVariable("WheelVoidThick",1.0);

  Control.addVariable("WheelInnerRadius",82.0);
  Control.addVariable("WheelCoolantRadius",125.0);
  Control.addVariable("WheelCaseRadius",127.0);
  Control.addVariable("WheelVoidRadius",130.0);
  Control.addVariable("WheelWMat",74001);
  Control.addVariable("WheelSteelMat",26317);
  Control.addVariable("WheelHeMat",2002);
  Control.addVariable("WheelInnerMat",26317);

  Control.addVariable("WheelNLayers",28);
  Control.addVariable("WheelRadius1",83.87);
  Control.addVariable("WheelMatTYPE1",1);
  Control.addVariable("WheelRadius2",85.25);
  Control.addVariable("WheelMatTYPE2",2);
  Control.addVariable("WheelRadius3",96.62);
  Control.addVariable("WheelMatTYPE3",3);
  Control.addVariable("WheelRadius4",97.12);
  Control.addVariable("WheelMatTYPE4",2);
  Control.addVariable("WheelRadius5",100.00);    // WRONG
  Control.addVariable("WheelMatTYPE5",3);
  Control.addVariable("WheelRadius6",101.97);
  Control.addVariable("WheelMatTYPE6",2);
  Control.addVariable("WheelRadius7",102.17);
  Control.addVariable("WheelMatTYPE7",1);
  Control.addVariable("WheelRadius8",102.57);
  Control.addVariable("WheelMatTYPE8",2);
  Control.addVariable("WheelRadius9",106.77);
  Control.addVariable("WheelMatTYPE9",3);
  Control.addVariable("WheelRadius10",107.09);
  Control.addVariable("WheelMatTYPE10",2);
  Control.addVariable("WheelRadius11",110.29);
  Control.addVariable("WheelMatTYPE11",3);
  Control.addVariable("WheelRadius12",110.58);
  Control.addVariable("WheelMatTYPE12",2); 
  Control.addVariable("WheelRadius13",112.58);
  Control.addVariable("WheelMatTYPE13",3); 
  Control.addVariable("WheelRadius14",113.3);
  Control.addVariable("WheelMatTYPE14",2); 
  Control.addVariable("WheelRadius15",114.45);
  Control.addVariable("WheelMatTYPE15",3); 
  Control.addVariable("WheelRadius16",114.7);
  Control.addVariable("WheelMatTYPE16",2); 
  Control.addVariable("WheelRadius17",114.9);
  Control.addVariable("WheelMatTYPE17",1); 
  Control.addVariable("WheelRadius18",116.5);
  Control.addVariable("WheelMatTYPE18",3); 
  Control.addVariable("WheelRadius19",116.75);
  Control.addVariable("WheelMatTYPE19",2); 
  Control.addVariable("WheelRadius20",118.35);
  Control.addVariable("WheelMatTYPE20",3); 
  Control.addVariable("WheelRadius21",118.6);
  Control.addVariable("WheelMatTYPE21",2); 
  Control.addVariable("WheelRadius22",120.0);
  Control.addVariable("WheelMatTYPE22",3); 
  Control.addVariable("WheelRadius23",120.25);
  Control.addVariable("WheelMatTYPE23",2); 
  Control.addVariable("WheelRadius24",121.65);
  Control.addVariable("WheelMatTYPE24",3); 
  Control.addVariable("WheelRadius25",121.9);
  Control.addVariable("WheelMatTYPE25",2); 
  Control.addVariable("WheelRadius26",123.1);
  Control.addVariable("WheelMatTYPE26",3); 
  Control.addVariable("WheelRadius27",123.35);
  Control.addVariable("WheelMatTYPE27",2); 
  Control.addVariable("WheelRadius28",124.55);
  Control.addVariable("WheelMatTYPE28",3); 


  // BE REFLECTOR
  Control.addVariable("BeRefXStep",0.0);  
  Control.addVariable("BeRefYStep",0.0);  
  Control.addVariable("BeRefZStep",0.0);
  Control.addVariable("BeRefXYangle",0.0); 
  Control.addVariable("BeRefZangle",0.0);
  Control.addVariable("BeRefRadius",30.0);
  Control.addVariable("BeRefHeight",90.0);
  Control.addVariable("BeRefWallThick",0.0);
  Control.addVariable("BeRefRefMat",4000);
  Control.addVariable("BeRefWallMat",0);
  Control.addVariable("BeRefFrontHeight",12.30);
 
 
  // FE REFLECTOR AND MONOLITH
  Control.addVariable("BulkXStep",0.0);
  Control.addVariable("BulkYStep",0.0);
  Control.addVariable("BulkZStep",0.0);
  Control.addVariable("BulkXYangle",0.0);
  Control.addVariable("BulkZangle",0.0);
  Control.addVariable("BulkNLayer",3);

  Control.Parse("BeRefRadius+0.0");
  Control.addVariable("BulkRadius1");

  //  Control.addVariable("BulkRadius1",47.5);
  Control.addVariable("BulkHeight1",45.0);
  Control.addVariable("BulkDepth1",45.0);
  Control.addVariable("BulkMat1",26316);

  Control.addVariable("BulkRadius2",65.0);
  Control.addVariable("BulkHeight2",75.0);
  Control.addVariable("BulkDepth2",75.0);
  Control.addVariable("BulkMat2",26316);           // stainless

  // Bulk steel layer [No individual guides]
  Control.addVariable("BulkRadius3",200.0);
  Control.addVariable("BulkHeight3",200.0);
  Control.addVariable("BulkDepth3",200.0);
  Control.addVariable("BulkMat3",26000);           // Bulk Steel


  // BULK FLIGHT VOID
  // Control.addVariable("BulkLAFlightSideIndex",-2);   // Index
  Control.addVariable("BulkLAFlightSideIndex",0);   // Index
  Control.addVariable("BulkLAFlightXStep",0.0);      // Step from centre
  Control.addVariable("BulkLAFlightZStep",0.0);      // Step from centre
  Control.addVariable("BulkLAFlightAngleXY1",30.0);  // Angle out
  Control.addVariable("BulkLAFlightAngleXY2",30.0);  // Angle out
  Control.addVariable("BulkLAFlightAngleZTop",0.0);  // Step down angle
  Control.addVariable("BulkLAFlightAngleZBase",0.0); // Step up angle
  Control.addVariable("BulkLAFlightHeight",10.0);    // Full height
  // Control.addVariable("BulkLAFlightWidth",23.3);     // Full width
  Control.addVariable("BulkLAFlightWidth",20.0);     // Full width
  Control.addVariable("BulkLAFlightNLiner",0);       // Liner


  // Guide BAY [ All 4 same ]
  Control.addVariable("GuideBayXStep",0.0);  
  Control.addVariable("GuideBayYStep",0.0);  
  Control.addVariable("GuideBayZStep",0.0);
  Control.addVariable("GuideBayZangle",0.0);
  Control.addVariable("GuideBayViewAngle",65.0); 
  Control.addVariable("GuideBayInnerHeight",11.0);
  Control.addVariable("GuideBayInnerDepth",11.0);
  Control.addVariable("GuideBayMidRadius",170.0);
  Control.addVariable("GuideBayHeight",32.0);
  Control.addVariable("GuideBayDepth",40.0);
  Control.addVariable("GuideBayMat",26000);
  // LL
  Control.addVariable("GuideBay1XYangle",0.0); 
  // LR
  Control.addVariable("GuideBay2XYangle",175.0); 
  // UL
  Control.addVariable("GuideBay3XYangle",2.5);
  // UR 
  Control.addVariable("GuideBay4XYangle",177.5); 
  Control.addVariable("GuideBay1NItems",12);  
  Control.addVariable("GuideBay2NItems",12);  
  Control.addVariable("GuideBay3NItems",12);  
  Control.addVariable("GuideBay4NItems",12);  


 // SHUTTER BAY
  Control.addVariable("ShutterBayXStep",0.0);  
  Control.addVariable("ShutterBayYStep",0.0);  
  Control.addVariable("ShutterBayZStep",0.0);
  Control.addVariable("ShutterBayXYangle",0.0); 
  Control.addVariable("ShutterBayZangle",0.0);
  Control.addVariable("ShutterBayRadius",600.0);
  Control.addVariable("ShutterBayHeight",400.0);
  Control.addVariable("ShutterBayDepth",400.0);
  Control.addVariable("ShutterBayMat",26000);


  EssBeamLinesVariables(Control);
  EssConicModerator(Control);
  EssWheel(Control);


  Control.addVariable("sdefEnergy",2500.0);  

  return;
}



void
EssConicModerator(FuncDataBase& Control)
  /*!
    Create all the Conic moderator option variables
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssConicModerator");
  
  // CONE MODERATOR
  Control.addVariable("LowConeModXStep",-0.0);      
  Control.addVariable("LowConeModYStep",4.0);      
  Control.addVariable("LowConeModZStep",-18.0);      
  Control.addVariable("LowConeModXYAngle",125.15);      
  Control.addVariable("LowConeModZAngle",0.0);      

  Control.addVariable("LowConeModIWidth",2.0);      
  Control.addVariable("LowConeModIHeight",2.0);      
  Control.addVariable("LowConeModOWidth",20.0);      
  Control.addVariable("LowConeModOHeight",10.0);      
  Control.addVariable("LowConeModLength",20.0);      
  Control.addVariable("LowConeModFaceThick",2.0);      
  Control.addVariable("LowConeModThick",1.5);      

  Control.addVariable("LowConeModAlThick",0.5);      

  Control.addVariable("LowConeModVacGap",0.3);      
  Control.addVariable("LowConeModWaterAlThick",0.5);      
  Control.addVariable("LowConeModWaterThick",2.0);      
  Control.addVariable("LowConeModVoidGap",0.3);      
  Control.addVariable("LowConeModWaterMat",11);      

  Control.addVariable("LowConeModModTemp",20.0);         // Temperature of H2 
  Control.addVariable("LowConeModModMat",25);            // Liquid H2
  Control.addVariable("LowConeModAlMat",5);              // Aluminium

  //  Control.addVariable("LowAFlightXStep",0.0);      // Step from centre


  // CONE MODERATOR
  Control.addVariable("LowConeModBXStep",-0.0);      
  Control.addVariable("LowConeModBYStep",4.0);      
  Control.addVariable("LowConeModBZStep",-18.0);      
  Control.addVariable("LowConeModBXYAngle",125.15-180.0);      
  Control.addVariable("LowConeModBZAngle",0.0);      

  Control.addVariable("LowConeModBIWidth",2.0);      
  Control.addVariable("LowConeModBIHeight",2.0);      
  Control.addVariable("LowConeModBOWidth",20.0);      
  Control.addVariable("LowConeModBOHeight",10.0);      
  Control.addVariable("LowConeModBLength",20.0);      
  Control.addVariable("LowConeModBFaceThick",2.0);      
  Control.addVariable("LowConeModBThick",1.5);      

  Control.addVariable("LowConeModBAlThick",0.5);      

  Control.addVariable("LowConeModBVacGap",0.3);      
  Control.addVariable("LowConeModBWaterAlThick",0.5);      
  Control.addVariable("LowConeModBWaterThick",2.0);      
  Control.addVariable("LowConeModBVoidGap",0.3);      
  Control.addVariable("LowConeModBWaterMat",11);      

  Control.addVariable("LowConeModBModTemp",20.0);         // Temperature of H2 
  Control.addVariable("LowConeModBModMat",25);            // Liquid H2
  Control.addVariable("LowConeModBAlMat",5);              // Aluminium

}


void
EssBeamLinesVariables(FuncDataBase& Control)
  /*!
    Create all the beamline variabes
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssBeamLinesVariables");
  for(int i=0;i<4;i++)
    {
      const std::string baseKey=
	StrFunc::makeString("G",i+1)+"BLine";
      // BeamLine in guide bay
      Control.addVariable(baseKey+"XStep",0.0);  
      Control.addVariable(baseKey+"YStep",0.0);  
      Control.addVariable(baseKey+"ZStep",0.0);
      Control.addVariable(baseKey+"Zangle",0.0);
      Control.addVariable(baseKey+"Mat",26316);
      Control.addVariable(baseKey+"BeamXYAngle",0.0); 
      Control.addVariable(baseKey+"BeamZAngle",0.0);
      Control.addVariable(baseKey+"BeamXStep",0.0);
      Control.addVariable(baseKey+"BeamZStep",0.0);
      Control.addVariable(baseKey+"BeamHeight",7.2);
      Control.addVariable(baseKey+"BeamWidth",7.6);
      Control.addVariable(baseKey+"NSegment",3);
      Control.addVariable(baseKey+"Width1",22.0);
      Control.addVariable(baseKey+"Height1",22.0);
      Control.addVariable(baseKey+"Width2",28.0);
      Control.addVariable(baseKey+"Height2",44.0);
      Control.addVariable(baseKey+"Width3",40.0);
      Control.addVariable(baseKey+"Height3",60.0);
      Control.addVariable(baseKey+"Length1",170.0);
      Control.addVariable(baseKey+"Length2",170.0);
      Control.addVariable(baseKey+"1XYangle",27.50); 
      Control.addVariable(baseKey+"2XYangle",22.5); 
      Control.addVariable(baseKey+"3XYangle",17.5); 
      Control.addVariable(baseKey+"4XYangle",12.5); 
      Control.addVariable(baseKey+"5XYangle",7.5); 
      Control.addVariable(baseKey+"6XYangle",2.5); 
      Control.addVariable(baseKey+"7XYangle",-2.5); 
      Control.addVariable(baseKey+"8XYangle",-7.5);
      Control.addVariable(baseKey+"9XYangle",-12.5); 
      Control.addVariable(baseKey+"10XYangle",-17.5); 
      Control.addVariable(baseKey+"11XYangle",-22.5); 
      Control.addVariable(baseKey+"12XYangle",-27.50); 
    }
  return;
}

}  // NAMESPACE setVariable
