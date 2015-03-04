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

#include <vector>

#include <opm/core/io/eclipse/EclipseWriteRFTHandler.hpp>
#include <opm/core/simulator/SimulatorState.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/miscUtilities.hpp>


#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellSet.hpp>

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>


#include <ert/ecl/ecl_rft_node.h>
#include <ert/ecl/ecl_rft_file.h>


namespace Opm {
namespace EclipseWriterDetails {

    EclipseWriteRFTHandler::EclipseWriteRFTHandler() {}

    void EclipseWriteRFTHandler::writeTimeStep(const std::string& filename,
                                               const SimulatorTimerInterface& simulatorTimer,
                                               std::vector<WellConstPtr>& wells,
                                               EclipseGridConstPtr eclipseGrid,
                                               size_t numCells,
                                               const int * compressedToCartesianCellIdx,
                                               std::vector<double>& pressure,
                                               std::vector<double>& swat,
                                               std::vector<double>& sgas) {



        std::vector<ecl_rft_node_type *> rft_nodes;
        for (std::vector<WellConstPtr>::const_iterator ci = wells.begin(); ci != wells.end(); ++ci) {
            WellConstPtr well = *ci;
            if ((well->getRFTActive(simulatorTimer.currentStepNum())) || (well->getPLTActive(simulatorTimer.currentStepNum()))) {
                ecl_rft_node_type * ecl_node = createEclRFTNode(well,
                                                                 simulatorTimer,
                                                                 eclipseGrid,
                                                                 numCells,
                                                                 compressedToCartesianCellIdx,
                                                                 pressure,
                                                                 swat,
                                                                 sgas);

                if (well->getPLTActive(simulatorTimer.currentStepNum())) {
                    std::cerr << "PLT not supported, writing RFT data" << std::endl;
                }

                rft_nodes.push_back(ecl_node);
            }
        }


        if (rft_nodes.size() > 0) {
            ecl_rft_file_update(filename.c_str(), rft_nodes.data(), rft_nodes.size(), ERT_ECL_METRIC_UNITS);
        }
    }




    ecl_rft_node_type * EclipseWriteRFTHandler::createEclRFTNode(WellConstPtr well,
                                                                  const SimulatorTimerInterface& simulatorTimer,
                                                                  EclipseGridConstPtr eclipseGrid,
                                                                  size_t numCells,
                                                                  const int * compressedToCartesianCellIdx,
                                                                  const std::vector<double>& pressure,
                                                                  const std::vector<double>& swat,
                                                                  const std::vector<double>& sgas) {


        const std::string& well_name      = well->name();
        size_t             timestep       = (size_t)simulatorTimer.currentStepNum();
        time_t             recording_date = simulatorTimer.currentPosixTime();
        double             days           = Opm::unit::convert::to(simulatorTimer.simulationTimeElapsed(), Opm::unit::day);

        std::vector<int> globalToActiveIndex = getGlobalToActiveIndex(compressedToCartesianCellIdx, numCells, eclipseGrid->getCartesianSize());

        std::string type = "RFT";
        ecl_rft_node_type * ecl_rft_node = ecl_rft_node_alloc_new(well_name.c_str(), type.c_str(), recording_date, days);

        CompletionSetConstPtr completionsSet = well->getCompletions(timestep);
        for (int index = 0; index < completionsSet->size(); ++index) {
            CompletionConstPtr completion = completionsSet->get(index);
            size_t i = (size_t)completion->getI();
            size_t j = (size_t)completion->getJ();
            size_t k = (size_t)completion->getK();

            size_t global_index = eclipseGrid->getGlobalIndex(i,j,k);
            int active_index = globalToActiveIndex[global_index];

            if (active_index > -1) {
                double depth = eclipseGrid->getCellDepth(i,j,k);
                double completion_pressure = pressure.size() > 0 ? pressure[active_index] : 0.0;
                double saturation_water    = swat.size() > 0 ? swat[active_index] : 0.0;
                double saturation_gas      = sgas.size() > 0 ? sgas[active_index] : 0.0;

                ecl_rft_cell_type * ecl_rft_cell = ecl_rft_cell_alloc_RFT( i ,j, k , depth, completion_pressure, saturation_water, saturation_gas);
                ecl_rft_node_append_cell( ecl_rft_node , ecl_rft_cell);
            }
        }

        return ecl_rft_node;
    }


    std::vector<int> EclipseWriteRFTHandler::getGlobalToActiveIndex(const int* compressedToCartesianCellIdx, size_t activeSize, size_t cartesianSize) {
        std::vector<int> globalToActiveIndex(cartesianSize, -1);
        for (int active_index = 0; active_index < activeSize; ++active_index) {
            int global_index = compressedToCartesianCellIdx[active_index];
            globalToActiveIndex[global_index] = active_index;
        }

        return globalToActiveIndex;
    }

}//namespace EclipseWriterDetails
}//namespace Opm
