/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   delftInc/RElement.h
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
#ifndef delftSystem_RElement_h
#define delftSystem_RElement_h

class FuncDataBase;

namespace delftSystem
{

/*!
  \class RElement
  \version 1.0
  \author S. Ansell
  \date June 2012
  \brief Element object of a reactors

  Holds a general element object of a reactor
*/

class RElement  : public attachSystem::FixedComp,
  public attachSystem::ContainedComp
{
 protected:

  const size_t XIndex;         ///< Index of block [X]
  const size_t YIndex;         ///< Index of block [Z]

  const int surfIndex;          ///< Index of surface offset
  int cellIndex;                ///< Cell index
 
  int insertCell;               ///< Cell to insert into

  double xStep;         ///< xStep
  double yStep;         ///< yStep
  double zStep;         ///< zStep
  double xyAngle;       ///< Rotation xy
  double zAngle;        ///< Rotation z

  void populate(const Simulation&);

 public:

  RElement(const size_t,const size_t,const std::string&);
  RElement(const RElement&);
  RElement& operator=(const RElement&);
  virtual ~RElement() {}   ///< Destructor

  virtual void createAll(Simulation&,const FixedComp&,
		 const Geometry::Vec3D&)=0;

};

}  

#endif
