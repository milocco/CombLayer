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
#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
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
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Triple.h"
#include "NRange.h"
#include "NList.h"
#include "Tally.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "surfDivide.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "KGroup.h"
#include "Source.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "generateSurf.h"
#include "SimProcess.h"
#include "chipDataStore.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "LinearComp.h"
#include "ContainedComp.h"
#include "channel.h"
#include "boxValues.h"
#include "boxUnit.h"
#include "BoxLine.h"


#include "ProtonTube.h"

namespace essSystem
{

ProtonTube::ProtonTube(const std::string& Key) :
  attachSystem::ContainedComp(),attachSystem::FixedComp(Key,3),
  ptIndex(ModelSupport::objectRegister::Instance().cell(Key)),
  cellIndex(ptIndex+1)

{}

ProtonTube::ProtonTube(const ProtonTube& A) :
  attachSystem::ContainedComp(A),attachSystem::FixedComp(A),
  ptIndex(A.ptIndex),cellIndex(A.cellIndex),
  xStep(A.xStep),yStep(A.yStep),zStep(A.zStep),
  xyAngle(A.xyAngle),zAngle(A.zAngle),nSecT(A.nSecT),
  radiusT(A.radiusT),lengthT(A.lengthT),cutT(A.cutT),
  thickT(A.thickT),inMatT(A.inMatT),wallMatT(A.wallMatT),
  XsetW(A.XsetW),YsetW(A.YsetW),ZsetW(A.ZsetW),sideW(A.sideW),
  nSecW(A.nSecW),radiusW(A.radiusW), thickW(A.thickW),matW(A.matW),
  frameSide(A.frameSide),frameHeightA(A.frameHeightA),
  frameHeightB(A.frameHeightB),frameHeightC(A.frameHeightC),
  frameHeightD(A.frameHeightD),  
  frameWidthA(A.frameWidthA),frameWidthB(A.frameWidthB),
  frameWidthC(A.frameWidthC),frameThickA(A.frameThickA),
  frameThickB(A.frameThickB),frameThickC(A.frameThickC),
  tubeN(A.tubeN),tubeRadius(A.tubeRadius),tubeThick(A.tubeThick),
  tubeHe(A.tubeHe),tubeAl(A.tubeAl),extTubeHe(A.extTubeHe),
  mainPBeamACell(A.mainPBeamACell)

  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}

ProtonTube&
ProtonTube::operator=(const ProtonTube& A) 
{
  if (&A!=this)
    {
     attachSystem::ContainedComp::operator=(A);
     attachSystem::FixedComp::operator=(A);
     //    wheelIndex=A.wheelIndex;
     cellIndex=A.cellIndex;
     xStep=A.xStep;
     yStep=A.yStep;
     zStep=A.zStep;
     xyAngle=A.xyAngle;
     zAngle=A.zAngle;
     nSecT=A.nSecT;
     radiusT=A.radiusT;
     lengthT=A.lengthT;
     cutT=A.cutT;
     thickT=A.thickT;
     inMatT=A.inMatT;
     wallMatT=A.wallMatT;
     XsetW=A.XsetW;
     YsetW=A.YsetW;
     ZsetW=A.ZsetW;
     sideW=A.sideW;
     nSecW=A.nSecW;
     radiusW=A.radiusW;
     thickW=A.thickW;
     matW=A.matW;    
     frameSide=A.frameSide;
     frameHeightA=A.frameHeightA;
     frameHeightB=A.frameHeightB;
     frameHeightC=A.frameHeightC;  
     frameHeightD=A.frameHeightD;  
     frameWidthA=A.frameWidthA;
     frameWidthB=A.frameWidthB;
     frameWidthC=A.frameWidthC;
     frameThickA=A.frameThickA;
     frameThickB=A.frameThickB;
     frameThickC=A.frameThickC;
     tubeN=A.tubeN;
     tubeRadius=A.tubeRadius;
     tubeThick=A.tubeThick;
     tubeHe=A.tubeHe;
     tubeAl=A.tubeAl;
     extTubeHe=A.extTubeHe;
     mainPBeamACell=A.mainPBeamACell;

    }
  return *this;
}


ProtonTube::~ProtonTube()
  /*!
    Destructor
   */
{}

void
ProtonTube::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("ProtonTube","populate");


    // Master values
  xStep=Control.EvalVar<double>(keyName+"XStep");
  yStep=Control.EvalVar<double>(keyName+"YStep");
  zStep=Control.EvalVar<double>(keyName+"ZStep");
  xyAngle=Control.EvalVar<double>(keyName+"XYangle");
  zAngle=Control.EvalVar<double>(keyName+"Zangle");
  nSecT=Control.EvalVar<size_t>(keyName+"SectionN");   
  double R(0.0);
  double L(0.0);
  double C(0.0);
  double T(0.0);
  int IM;
  int WM;

  for(size_t i=0;i<nSecT;i++)
    {
     R=Control.EvalVar<double>
	(StrFunc::makeString(keyName+"Radius",i+1));   
     L+=Control.EvalVar<double>
	(StrFunc::makeString(keyName+"Length",i+1));   
     C=Control.EvalVar<double>
	(StrFunc::makeString(keyName+"Zcut",i+1));   
     T=Control.EvalVar<double>
     (StrFunc::makeString(keyName+"WallThick",i+1));     
     IM=Control.EvalVar<int>
	(StrFunc::makeString(keyName+"InnerMat",i+1));   
     WM=Control.EvalVar<int>
	(StrFunc::makeString(keyName+"WallMat",i+1));   

      radiusT.push_back(R);
      lengthT.push_back(L);
      cutT.push_back(C);
      thickT.push_back(T);
      inMatT.push_back(IM);
      wallMatT.push_back(WM);
    }

  XsetW=Control.EvalVar<double>(keyName+"WindowOffsetX");   
  YsetW=Control.EvalVar<double>(keyName+"WindowOffsetY");   
  ZsetW=Control.EvalVar<double>(keyName+"WindowOffsetZ");   
  sideW=Control.EvalVar<double>(keyName+"WindowBoxSide");   

  nSecW=Control.EvalVar<size_t>(keyName+"WindowBoxSectionsN");   
  double RW(0.0);
  double TW(0.0);
  int MW;

  for(size_t i=0;i<nSecW;i++)
    {
     RW=Control.EvalVar<double>
	(StrFunc::makeString(keyName+"WindowBoxRadius",i+1));   
     TW+=Control.EvalVar<double>
	(StrFunc::makeString(keyName+"WindowBoxThick",i+1));     
     MW=Control.EvalVar<int>
	(StrFunc::makeString(keyName+"WindowBoxMat",i+1));   

      radiusW.push_back(RW);
      thickW.push_back(TW);
      matW.push_back(MW);
    }
 
  frameSide=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxSide5"));    
  frameHeightA=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxHeightA5"));    
  frameHeightB=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxHeightB5"));    
  frameHeightC=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxHeightC5"));    
  frameHeightD=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxHeightD5"));    
  frameWidthA=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxWidthA5"));    
  frameWidthB=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxWidthB5"));    
  frameWidthC=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxWidthC5"));    
  frameThickA=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxThickA5"));    
  frameThickB=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxThickB5"));    
  frameThickC=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxThickC5"));    

  tubeN=Control.EvalVar<int> (StrFunc::makeString(keyName+"WindowBoxTubeN"));      
  tubeRadius=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxTubeRadius"));    
  tubeThick=Control.EvalVar<double>(StrFunc::makeString(keyName+"WindowBoxTubeThick"));    
  
  tubeHe=Control.EvalVar<int>(keyName+"WindowBoxTubeHeMat");   
  tubeAl=Control.EvalVar<int>(keyName+"WindowBoxTubeAlMat");  
  extTubeHe=Control.EvalVar<int>(keyName+"WindowBoxExtHeMat");

  return;
}

void
ProtonTube::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: Fixed Component
  */
{
  ELog::RegMethod RegA("ProtonTube","createUnitVector");
  attachSystem::FixedComp::createUnitVector(FC);
  
  
  Origin+=X*xStep+Y*(yStep-115)+Z*zStep;
  const Geometry::Quaternion Qz=
    Geometry::Quaternion::calcQRotDeg(zAngle,X);
  const Geometry::Quaternion Qxy=
    Geometry::Quaternion::calcQRotDeg(xyAngle,Z);
  Qz.rotate(Y);
  Qz.rotate(Z);
  Qxy.rotate(Y);
  Qxy.rotate(X);
  Qxy.rotate(Z);

  return;
}

void
ProtonTube::createSurfaces()
  /*!

  */
{
  ELog::RegMethod RegA("ProtonTube","createSurfaces");
  
  ModelSupport::buildCylinder(SMap,ptIndex+700,Origin,Z,600);  
  ModelSupport::buildPlane(SMap,ptIndex+600,Origin,Y);  
    


  int PT(ptIndex);
  for(size_t i=0;i<nSecT;i++)
    {
     ModelSupport::buildCylinder(SMap,PT+7,Origin,Y,radiusT[i]);  
     ModelSupport::buildCylinder(SMap,PT+8,Origin,Y,radiusT[i]+thickT[i]); 
     ModelSupport::buildPlane(SMap,PT+3,Origin-Y*lengthT[i],Y);  
     ModelSupport::buildPlane(SMap,PT+5,Origin-Z*(radiusT[i]-cutT[i]),Z);
     ModelSupport::buildPlane(SMap,PT+6,Origin+Z*(radiusT[i]-cutT[i]),Z);  
     PT+=10;
    }

  int PW(ptIndex+1000);
  ModelSupport::buildPlane(SMap,PW+1,Origin+X*(XsetW-sideW/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,PW+2,Origin+X*(XsetW+sideW/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,PW+3,Origin+X*XsetW+Y*(YsetW+sideW/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,PW+4,Origin+X*XsetW+Y*(YsetW-sideW/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,PW+5,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-sideW/2.0),Z);  
  ModelSupport::buildPlane(SMap,PW+6,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+sideW/2.0),Z);  

  for(size_t i=0;i<nSecW+3;i++)
    {
     PW+=10;
     if(i>=nSecW)
       { 
         radiusW[i]=radiusW[nSecW-(i-nSecW+2)];
         thickW[i]=thickW[i-1]+(thickW[nSecW-(i-nSecW+2)]-thickW[nSecW-(i-nSecW+3)]);
         matW[i]=matW[nSecW-(i-nSecW+2)];
	   }
     if(i==nSecW+2) thickW[i]+=thickW[0];
     ModelSupport::buildCylinder(SMap,PW+7,Origin,Y,radiusW[i]);  
     ModelSupport::buildPlane(SMap,PW+3,Origin+Y*(YsetW+sideW/2.0-thickW[i]),Y);  
    }

  int TI(ptIndex+2000);

  ModelSupport::buildPlane(SMap,TI+1,Origin+X*(XsetW-frameSide/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+2,Origin+X*(XsetW+frameSide/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+3,Origin+X*XsetW+Y*(YsetW-frameSide/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+4,Origin+X*XsetW+Y*(YsetW+frameSide/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+5,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-frameSide/2.0),Z);  
  ModelSupport::buildPlane(SMap,TI+6,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+frameSide/2.0),Z);  

  ModelSupport::buildPlane(SMap,TI+11,Origin+X*(XsetW-frameWidthA/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+12,Origin+X*(XsetW+frameWidthA/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+13,Origin+X*XsetW+Y*(YsetW-frameThickA/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+14,Origin+X*XsetW+Y*(YsetW+frameThickA/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+15,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-frameHeightA/2.0),Z);  
  ModelSupport::buildPlane(SMap,TI+16,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+frameHeightA/2.0),Z);  

  ModelSupport::buildPlane(SMap,TI+21,Origin+X*(XsetW-frameWidthB/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+22,Origin+X*(XsetW+frameWidthB/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+23,Origin+X*XsetW+Y*(YsetW-frameThickB/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+24,Origin+X*XsetW+Y*(YsetW+frameThickB/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+25,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-frameHeightB/2.0),Z);  
  ModelSupport::buildPlane(SMap,TI+26,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+frameHeightB/2.0),Z);  

  ModelSupport::buildPlane(SMap,TI+31,Origin+X*(XsetW-frameWidthC/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+32,Origin+X*(XsetW+frameWidthC/2.0)+Y*YsetW+Z*ZsetW,X);  
  ModelSupport::buildPlane(SMap,TI+33,Origin+X*XsetW+Y*(YsetW-frameThickC/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+34,Origin+X*XsetW+Y*(YsetW+frameThickC/2.0)+Z*ZsetW,Y);  
  ModelSupport::buildPlane(SMap,TI+35,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-frameHeightC/2.0),Z);  
  ModelSupport::buildPlane(SMap,TI+36,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+frameHeightC/2.0),Z);  

  ModelSupport::buildPlane(SMap,TI+45,Origin+X*XsetW+Y*YsetW+Z*(ZsetW-frameHeightD/2.0),Z);  
  ModelSupport::buildPlane(SMap,TI+46,Origin+X*XsetW+Y*YsetW+Z*(ZsetW+frameHeightD/2.0),Z);  


  for(int i=0;i<tubeN;i++)
    {
    ModelSupport::buildCylinder(SMap,TI+7,
                  Origin+X*(XsetW-frameWidthC/2.0+0.3+i*2*(tubeRadius+tubeThick))
                 +Y*YsetW+Z*ZsetW,Z,tubeRadius);  
    ModelSupport::buildCylinder(SMap,TI+8,
                  Origin+X*(XsetW-frameWidthC/2.0+0.3+i*2*(tubeRadius+tubeThick))
                  +Y*YsetW+Z*ZsetW,Z,tubeRadius+tubeThick);
    TI+=10;

    }
  return; 
}

void
ProtonTube::addToInsertChain(attachSystem::ContainedComp& CC) const
  /*!
    Adds this object to the containedComp to be inserted.
    \param CC :: ContainedComp object to add to this
  */
{
  for(int i=ptIndex+1;i<cellIndex;i++)
    CC.addInsertCell(i);
    
  return;
}

void
ProtonTube::createObjects(Simulation& System, const std::string& TargetSurfBoundary)
  /*!
    Adds the components
    \param System :: Simulation to create objects in
    \param TargetSurfBoundary :: boundary layer [expect to be target edge]
    \param RefSurfBoundary :: boundary layer [expect to be reflector edge]
  */
{
  ELog::RegMethod RegA("ProtonVoid","createObjects");
    int protonVoidCell;           ///< Inner void cell

  std::string OutTot;
  std::string OutTube;

  std::string OutLeft;
  OutLeft=TargetSurfBoundary;
  std::string Out;
  int PT(ptIndex);
  for(size_t i=0;i<nSecT;i++)
    {
      if (i==0) 
      {
       Out=ModelSupport::getComposite(SMap,ptIndex,PT, "(-7M 3M 5 -6 -600)");
       Out+=TargetSurfBoundary;
       System.addCell(MonteCarlo::Qhull(cellIndex++,inMatT[i],0.0,Out));
       OutTot=Out;
       }
      else
      {
        Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(-7M -3M 13M)");
        if(i==nSecT-1) Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(-1053 -7M -700)");
        System.addCell(MonteCarlo::Qhull(cellIndex++,inMatT[i],0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(7M -8M -3M 13M)");
        if(i==nSecT-1) Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(-1053 7M -8M -700)");
        System.addCell(MonteCarlo::Qhull(cellIndex++,wallMatT[i],0.0,Out));
        Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(-8M -3M 13M)");
        if(i==nSecT-1) Out=ModelSupport::getComposite(SMap,ptIndex,PT,"(-1053 -8M -700)");
       OutTot+=":"+Out;
       PT+=10; 
      }
    }

  int PW(ptIndex+1000);
  for(size_t i=0;i<nSecW+3;i++)
    {
     Out=ModelSupport::getComposite(SMap,ptIndex,PW,"(38 -17M -3M 13M)");


     if(i==nSecW-1) 
       {
        Out=ModelSupport::getComposite(SMap,ptIndex,"( 38 ");
        Out+=ModelSupport::getComposite(SMap,ptIndex+1000,PW," 1 -2 -3 13M 5 -6 )");
        OutTot+=":"+Out;  
        Out=ModelSupport::getComposite(SMap,ptIndex,"( 38 ");
        Out+=ModelSupport::getComposite(SMap,ptIndex+1000,PW," 1 -2 4 -3M 5 -6 )");
        OutTot+=":"+Out;  
        Out=ModelSupport::getComposite(SMap,ptIndex+1000,PW,"( 1 -2 -3M 13M 5 -6 )");
        OutTot+=":"+Out;  
	// 
        int TI(ptIndex+2000);

        for(int i=0;i<tubeN;i++)
          {
	  Out=ModelSupport::getComposite(SMap,ptIndex+2000,TI,"-7M 25 -26");
         System.addCell(MonteCarlo::Qhull(cellIndex++,tubeHe,0.0,Out));
         Out=ModelSupport::getComposite(SMap,ptIndex+2000,TI,"7M -8M 25 -26");
         System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out));

         OutTube+=ModelSupport::getComposite(SMap,TI," 8M ");

         TI+=10;  
           }
	// hugly piece by piece

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW,
                                       "(-11:12:-15:16) 1 -2 5 -6 -3M 13M");
        System.addCell(MonteCarlo::Qhull(cellIndex++,matW[i],0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 11 -12 13 -14 15 -16 (-25:26)");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeHe,0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW," 11 -12 -3M 14 15 -16 (-35:36) ");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW," 11 -12 13M -13 15 -16 (-35:36) ");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000," 11 -12 13 -14 25 -26 (-45:46)");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out+OutTube));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000," (11 -12 33 -34  45 -46) ");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out+OutTube));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 31 -32 -33 13 45 -46");
        System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out+OutTube));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 11 -12 (-31: 32) -33 13 45 -46");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out+OutTube));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 31 -32 34 -14 45 -46");
        System.addCell(MonteCarlo::Qhull(cellIndex++,extTubeHe,0.0,Out+OutTube));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW,"11 -12 (-31 : 32) 34 -14 45 -46");
        System.addCell(MonteCarlo::Qhull(cellIndex++,tubeAl,0.0,Out+OutTube));

       Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 11 -12 13M -13 35 -36");
       System.addCell(MonteCarlo::Qhull(cellIndex++,0,0.0,Out));

        Out=ModelSupport::getComposite(SMap,ptIndex+2000,PW, " 11 -12 -3M 14 35 -36");
        System.addCell(MonteCarlo::Qhull(cellIndex++,extTubeHe,0.0,Out));
      }
     else
       {
        System.addCell(MonteCarlo::Qhull(cellIndex++,matW[i],0.0,Out));
       }

     Out=ModelSupport::getComposite(SMap,ptIndex+1000," 1 -2 5 -6 ");

     if(i==nSecW-1) 
       {
	 Out+=ModelSupport::getComposite(SMap,ptIndex+2000,PW,"(-1:2:-5:6)-3M 13M ");
       }
     else
       {
     Out+=ModelSupport::getComposite(SMap,PW,"17M -3M 13M ");
       }
     System.addCell(MonteCarlo::Qhull(cellIndex++,inMatT[0],0.0,Out));

     PW+=10; 
     }
    addOuterSurf(OutTot);
 
  return;
}


void
ProtonTube::createLinks()
  /*!
    Creates a full attachment set
  */
{
  
  return;
}


void
ProtonTube::createAll(Simulation& System,
		      const attachSystem::FixedComp& TargetFC,
		      const long int tIndex)
  /*!
    Global creation of the hutch
    \param System :: Simulation to add vessel to
    \param FC :: FixedComp for origin
    \param tIndex :: Target plate surface
    for me target surf is enough
  */
{
  ELog::RegMethod RegA("ProtonTube","createAll");
  //  populate(System);
  populate(System.getDataBase());

  createUnitVector(TargetFC);
  createSurfaces();
  const std::string TSurf=(tIndex<0) ? 
    TargetFC.getLinkComplement(static_cast<size_t>(-(tIndex+1))) : 
    TargetFC.getLinkString(static_cast<size_t>(tIndex));

    createObjects(System,TSurf);
    createLinks();
    insertObjects(System); 
  return;
}

}  // NAMESPACE instrumentSystem
