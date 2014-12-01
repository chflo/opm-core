
/*
  Copyright 2014 Statoil IT
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
#include "config.h"

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE EclipseWriter
#include <boost/test/unit_test.hpp>

#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

// ERT stuff
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl_well/well_info.h>
#include <ert/ecl_well/well_state.h>
#include <ert/util/test_work_area.h>

#include <string.h>


void verifyWellState(const std::string& rst_filename,
                     Opm::EclipseGridConstPtr ecl_grid,
                     Opm::ScheduleConstPtr schedule) {

  well_info_type * well_info = well_info_alloc(ecl_grid->c_ptr());
  well_info_load_rstfile(well_info, rst_filename.c_str(), false);

  int numwells = well_info_get_num_wells(well_info);
  BOOST_CHECK(numwells == (int)schedule->numWells());

  std::vector<Opm::WellConstPtr> wells = schedule->getWells();

  for (int i = 0; i < numwells; ++i) {

    //Verify wellnames
    const char * wellname = well_info_iget_well_name(well_info, i);
    Opm::WellConstPtr well = wells.at(i);
    BOOST_CHECK(wellname == well->name());

    // Verify well-head position data
    well_ts_type * well_ts = well_info_get_ts(well_info , wellname);
    well_state_type * well_state = well_ts_iget_state(well_ts, 0);
    const well_conn_type * well_head = well_state_get_wellhead(well_state, ECL_GRID_GLOBAL_GRID);
    BOOST_CHECK(well_conn_get_i(well_head) == well->getHeadI());
    BOOST_CHECK(well_conn_get_j(well_head) == well->getHeadJ());

    //Verify number of completion connections
    for (int i = 0; i < well_ts_get_size(well_ts); ++i) {
      well_state_type * well_state = well_ts_iget_state(well_ts, i);
      const well_conn_collection_type * well_connections = well_state_get_global_connections( well_state );
      int num_wellconnections = well_conn_collection_get_size(well_connections);

      int report_nr = well_state_get_report_nr(well_state);
      Opm::CompletionSetConstPtr completions_set = well->getCompletions((size_t)report_nr);

      BOOST_CHECK(num_wellconnections == completions_set->size());
    }
  }

  well_info_free(well_info);
}


std::shared_ptr<Opm::BlackoilState> createBlackOilState(Opm::EclipseGridConstPtr eclGrid) {

  std::shared_ptr<Opm::GridManager> ourFineGridManagerPtr(new Opm::GridManager(eclGrid));
  std::shared_ptr<Opm::BlackoilState> blackoilState(new Opm::BlackoilState);
  blackoilState->init(*ourFineGridManagerPtr->c_grid(), 3);

  return blackoilState;
}


Opm::DeckConstPtr createDeck(const std::string& eclipse_data_filename) {
  Opm::ParserPtr parser(new Opm::Parser());
  Opm::ParserLogPtr parserLog(new Opm::ParserLog);
  Opm::DeckConstPtr deck = parser->parseFile(eclipse_data_filename, true, parserLog);

  return deck;
}


Opm::EclipseWriterPtr createEclipseWriter(Opm::DeckConstPtr deck,
                                          Opm::EclipseStatePtr eclipseState,
                                          std::string& eclipse_data_filename) {

  Opm::parameter::ParameterGroup params;
  params.insertParameter("deck_filename", eclipse_data_filename);

  const Opm::PhaseUsage phaseUsage = Opm::phaseUsageFromDeck(deck);

  Opm::EclipseWriterPtr eclWriter(new Opm::EclipseWriter(params,
                                                         eclipseState,
                                                         phaseUsage,
                                                         eclipseState->getEclipseGrid()->getCartesianSize(),
                                                         0));
  return eclWriter;
}


BOOST_AUTO_TEST_CASE(EclipseWriteRestartWellInfo)
{
    std::string eclipse_data_filename    = "testBlackoilState3.DATA";
    std::string eclipse_restart_filename = "TESTBLACKOILSTATE3.UNRST";

    test_work_area_type * test_area = test_work_area_alloc("TEST_EclipseWriteNumWells");
    test_work_area_copy_file(test_area, eclipse_data_filename.c_str());

    Opm::DeckConstPtr     deck = createDeck(eclipse_data_filename);
    Opm::EclipseStatePtr  eclipseState(new Opm::EclipseState(deck));
    Opm::EclipseWriterPtr eclipseWriter = createEclipseWriter(deck, eclipseState, eclipse_data_filename);

    std::shared_ptr<Opm::SimulatorTimer> simTimer( new Opm::SimulatorTimer() );
    simTimer->init(eclipseState->getSchedule()->getTimeMap());

    eclipseWriter->writeInit(*simTimer);

    std::shared_ptr<Opm::WellState> wellState(new Opm::WellState());
    std::shared_ptr<Opm::BlackoilState> blackoilState = createBlackOilState(eclipseState->getEclipseGrid());
    wellState->init(0, *blackoilState);

    int countTimeStep = eclipseState->getSchedule()->getTimeMap()->numTimesteps();

    for(int timestep=0; timestep <= countTimeStep; ++timestep){
      simTimer->setCurrentStepNum(timestep);
      eclipseWriter->writeTimeStep(*simTimer, *blackoilState, *wellState);
    }

    verifyWellState(eclipse_restart_filename, eclipseState->getEclipseGrid(), eclipseState->getSchedule());

    test_work_area_free(test_area);
}
