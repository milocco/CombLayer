/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   t1UpgradeInc/ConicInfo.h
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
#ifndef ts1System_ConicInfo_h
#define ts1System_ConicInfo_h

class Simulation;

namespace ts1System
{

  /*!
    \class ConicInfo
   */
class ConicInfo
{
 private:
 
  int cylFlag;              ///< Cylinder type / Cone type
  Geometry::Vec3D Cent;     ///< Conic
  Geometry::Vec3D Axis;     ///< Axis
  double angle;             ///< Angle [deg] / Radius
  double wall;              ///< Wall thickness
  int wallMat;              ///< Material

 public:

  ConicInfo(const Geometry::Vec3D&,const Geometry::Vec3D&,
	    const double,const double,const int);
  ConicInfo(const Geometry::Vec3D&,const Geometry::Vec3D&,
	    const double,const double,const int,const int);
  ConicInfo(const ConicInfo&);
  ConicInfo& operator=(const ConicInfo&);
  ~ConicInfo() {}  ///< Destructor

  Geometry::Vec3D getCent(const Geometry::Vec3D&,
			  const Geometry::Vec3D&,
			  const Geometry::Vec3D&) const;
  Geometry::Vec3D getAxis(const Geometry::Vec3D&,
			  const Geometry::Vec3D&,
			  const Geometry::Vec3D&) const;

  /// Is cylinder
  int isCyl() const { return cylFlag; }
  /// Wall material
  int getMat() const { return wallMat; } 
  /// Wall thicknes
  double getWall() const { return wall; }  
  /// Angle
  double getAngle() const { return angle; }  
  /// get Radius
  double getRadius() const { return angle; }  
};

}

#endif
 
