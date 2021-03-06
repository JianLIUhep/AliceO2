# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(
  MathUtils
  SOURCES src/CachingTF1.cxx
          src/Cartesian2D.cxx
          src/Cartesian3D.cxx
          src/Chebyshev3D.cxx
          src/Chebyshev3DCalc.cxx
          src/MathBase.cxx
          src/RandomRing.cxx
          src/Primitive2D.cxx
  PUBLIC_LINK_LIBRARIES
    ROOT::Hist
    FairRoot::Base
    O2::CommonConstants
    O2::GPUCommon
    ROOT::GenVector
    ROOT::Geom
    Vc::Vc)

o2_target_root_dictionary(
  MathUtils
  HEADERS include/MathUtils/Utils.h
          include/MathUtils/Chebyshev3D.h
          include/MathUtils/Chebyshev3DCalc.h
          include/MathUtils/MathBase.h
          include/MathUtils/Cartesian2D.h
          include/MathUtils/Cartesian3D.h
          include/MathUtils/CachingTF1.h
          include/MathUtils/RandomRing.h
          include/MathUtils/Primitive2D.h
          include/MathUtils/Bracket.h)

o2_add_test(
  CachingTF1
  SOURCES test/testCachingTF1.cxx
  COMPONENT_NAME MathUtils
  PUBLIC_LINK_LIBRARIES O2::MathUtils
  LABELS utils)

o2_add_test(
  Cartesian3D
  SOURCES test/testCartesian3D.cxx
  COMPONENT_NAME MathUtils
  PUBLIC_LINK_LIBRARIES O2::MathUtils
  LABELS utils)

o2_add_test(
  Utils
  SOURCES test/testUtils.cxx
  COMPONENT_NAME MathUtils
  PUBLIC_LINK_LIBRARIES O2::MathUtils
  LABELS utils)
