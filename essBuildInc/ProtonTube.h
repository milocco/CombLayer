#ifndef essSystem_ProtonTube_h
#define essSystem_ProtonTube_h

class Simulation;

namespace essSystem
{

/*!
  \class ProtonTube
  \author A. Milocco
  \version 1.0
  \date February 2013
  \brief Reflector object 
*/

class ProtonTube : public attachSystem::ContainedComp,
    public attachSystem::FixedComp
{
 private:

  const int ptIndex;             ///< Index of surface offset
  int cellIndex;                  ///< Cell index

  double xStep;                   ///< X step
  double yStep;                   ///< Y step
  double zStep;                   ///< Z step
  double xyAngle;                 ///< XY Angle
  double zAngle;                  ///< Z Angle
  size_t nSecT;
  std::vector<double> radiusT;    
  std::vector<double> lengthT;        
  std::vector<double> cutT;        
  std::vector<double> thickT;        
  std::vector<int> inMatT;     
  std::vector<int> wallMatT;     

  double XsetW;                 
  double YsetW;                 
  double ZsetW;                 
  double sideW;                 
  size_t nSecW;
  std::vector<double> radiusW;    
  std::vector<double> thickW;        
  std::vector<int> matW;     
  double frameSide;
  double frameHeightA;
  double frameHeightB;
  double frameHeightC;
  double frameHeightD;
  double frameWidthA;
  double frameWidthB;
  double frameWidthC;
  double frameThickA;  
  double frameThickB;
  double frameThickC;
  int tubeN;
  double tubeRadius;
  double tubeThick;
  int tubeHe;
  int tubeAl;
  int extTubeHe;

  int mainPBeamACell; 


  // Functions:

  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&);
  void createSurfaces();
  void createObjects(Simulation&,const std::string&);
  void createLinks();

 public:

  ProtonTube(const std::string&);
  ProtonTube(const ProtonTube&);
  ProtonTube& operator=(const ProtonTube&);
  virtual ~ProtonTube();
   

  int getCell() const { return mainPBeamACell; }
  virtual void addToInsertChain(attachSystem::ContainedComp&) const; 
   void createAll(Simulation&,const attachSystem::FixedComp&, const long int);
 
};

}

#endif
 
