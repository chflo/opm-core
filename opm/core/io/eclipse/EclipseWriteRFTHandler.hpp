/*
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP
#define OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP

#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorState.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <ert/ecl/ecl_rft_node.h>


namespace Opm {
namespace EclipseWriterDetails {


    class EclipseWriteRFTHandler {

    public:
    EclipseWriteRFTHandler();


    static void writeTimeStep(const std::string& filename,
                              const SimulatorTimerInterface& simulatorTimer,
                              std::vector<WellConstPtr>& wells,
                              EclipseGridConstPtr eclipseGrid,
                              size_t numCells,
                              const int * compressedToCartesianCellIdx,
                              std::vector<double>& pressure,
                              std::vector<double>& swat,
                              std::vector<double>& sgas);



    private:

    static ecl_rft_node_type * createEclRFTNode(WellConstPtr well,
                                                 const SimulatorTimerInterface& simulatorTimer,
                                                 EclipseGridConstPtr eclipseGrid,
                                                 size_t numCells,
                                                 const int * compressedToCartesianCellIdx,
                                                 const std::vector<double>& pressure,
                                                 const std::vector<double>& swat,
                                                 const std::vector<double>& sgas);

    static std::vector<int> getGlobalToActiveIndex(const int * compressedToCartesianCellIdx, size_t activeSize, size_t cartesianSize);

    };




}//namespace EclipseWriterDetails
}//namespace Opm


#endif //OPM_ECLIPSE_WRITE_RFT_HANDLER_HPP
