
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
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

// ERT stuff
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/util/test_work_area.h>

#include <string.h>


void verifyNumWells(const std::string&  filename, int timestep, int numWellsFromSchedule){
  ecl_file_type * restart_file = ecl_file_open(filename.c_str(), 0);
  ecl_file_load_all(restart_file);

  const ecl_kw_type * intehead_kw = ecl_file_iget_named_kw(restart_file, INTEHEAD_KW, timestep);
  int numwells = ecl_kw_iget_int(intehead_kw, INTEHEAD_NWELLS_INDEX);
  BOOST_ASSERT(numWellsFromSchedule == numwells);

  ecl_file_close(restart_file);
}


void verifyIwelData(const std::string& filename, int timestep) {
  ecl_file_type * restart_file = ecl_file_open(filename.c_str(), 0);
  ecl_file_load_all(restart_file);

  const ecl_kw_type * intehead_kw = ecl_file_iget_named_kw(restart_file, INTEHEAD_KW, timestep);
  int numwells = ecl_kw_iget_int(intehead_kw, INTEHEAD_NWELLS_INDEX);
  int NIWELZ = ecl_kw_iget_int(intehead_kw, INTEHEAD_NIWELZ_INDEX);
  int expected_IWEL_data_size = numwells * NIWELZ;

  const ecl_kw_type * iwel_kw = ecl_file_iget_named_kw(restart_file, IWEL_KW, timestep);
  BOOST_ASSERT(expected_IWEL_data_size == ecl_kw_get_size(iwel_kw));

  ecl_file_close(restart_file);
}


void verifyZwelData(const std::string&  filename, int timestep) {
  ecl_file_type * restart_file = ecl_file_open(filename.c_str(), 0);
  ecl_file_load_all(restart_file);

  const ecl_kw_type * intehead_kw = ecl_file_iget_named_kw(restart_file, INTEHEAD_KW, timestep);
  int numwells = ecl_kw_iget_int(intehead_kw, INTEHEAD_NWELLS_INDEX);
  int NZWELZ = ecl_kw_iget_int(intehead_kw, INTEHEAD_NZWELZ_INDEX);
  int expected_ZWEL_data_size = numwells * NZWELZ;

  const ecl_kw_type * zwel_kw = ecl_file_iget_named_kw(restart_file, ZWEL_KW, timestep);
  BOOST_ASSERT(expected_ZWEL_data_size == ecl_kw_get_size(zwel_kw));

  ecl_file_close(restart_file);
}


void  verifyIconData(const std::string&  filename, int timestep) {
  ecl_file_type * restart_file = ecl_file_open(filename.c_str(), 0);
  ecl_file_load_all(restart_file);

  const ecl_kw_type * intehead_kw = ecl_file_iget_named_kw(restart_file, INTEHEAD_KW, timestep);
  int numwells = ecl_kw_iget_int(intehead_kw, INTEHEAD_NWELLS_INDEX); //Num wells
  int NICONZ   = ecl_kw_iget_int(intehead_kw, INTEHEAD_NICONZ_INDEX); //Num elements per connection element
  int NCWMAX   = ecl_kw_iget_int(intehead_kw, INTEHEAD_NCWMAX_INDEX); //Max number of connections for a well

  int expected_ICON_data_size = NICONZ * NCWMAX * numwells;

  const ecl_kw_type * icon_kw = ecl_file_iget_named_kw(restart_file, ICON_KW, timestep);
  BOOST_ASSERT(expected_ICON_data_size == ecl_kw_get_size(icon_kw));

  ecl_file_close(restart_file);
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

    for (int timestep = 0; timestep <= countTimeStep; ++timestep) {
      int numWellsFromSchedule  = eclipseState->getSchedule()->numWells(timestep);
      verifyNumWells(eclipse_restart_filename, timestep, numWellsFromSchedule);
      verifyIwelData(eclipse_restart_filename, timestep);
      verifyZwelData(eclipse_restart_filename, timestep);
      verifyIconData(eclipse_restart_filename, timestep);
    }

    test_work_area_free(test_area);
}
