/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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


#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif

#define BOOST_TEST_MODULE WellStateResizeTest

#include <opm/core/wells.h>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/core/simulator/BlackoilState.hpp>

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <memory>

using namespace Opm;
/*
Example of COORD and ZCORN keywords found in opm-parser/opm/parser/eclipse/EclipseState/Grid/tests/EclipseGridTests.cpp.
input2 (with DXV, DYV, DZV and TOPS keywords) found in opm-core/tests/test_EclipseWriteRFTHandler.cpp.
If problem with the input data deck, the file opm-core/tests/wells_manager_data.data seems to have data that works
(gives us wells in c_wells() in the wellsManager).
*/

void testSizes(WellState& wsResize, size_t xpectedNumWells, size_t xpectedTotalNumCompl, size_t numPhases){
    size_t szBhp = wsResize.bhp().size();
    size_t szTemperature = wsResize.temperature().size();
    size_t szWellRates = wsResize.wellRates().size();
    size_t szPerfRates = wsResize.perfRates().size();
    size_t szPerfPress = wsResize.perfPress().size();

    BOOST_CHECK_EQUAL(xpectedNumWells, szBhp);
    BOOST_CHECK_EQUAL(xpectedNumWells, szTemperature);
    BOOST_CHECK_EQUAL((xpectedNumWells * numPhases), szWellRates);
    BOOST_CHECK_EQUAL(xpectedTotalNumCompl, szPerfRates);
    BOOST_CHECK_EQUAL(xpectedTotalNumCompl, szPerfPress);
}

BOOST_AUTO_TEST_CASE(resizeWellState) {
    Opm::Parser parser;

    std::string input =
            "RUNSPEC\n"
            "OIL\n"
            "GAS\n"
            "WATER\n"
            "DIMENS\n"
            " 10 10 10 /\n"

            "GRID\n"
            "DXV\n"
            "10*0.25 /\n"
            "DYV\n"
            "10*0.25 /\n"
            "DZV\n"
            "10*0.25 /\n"
            "TOPS\n"
            "100*0.25 /\n"
            "\n"

            "START             -- 0 \n"
            "1 NOV 1979 / \n"

            "SCHEDULE\n"
            "DATES             -- 1\n"
            " 10  OKT 2008 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_1'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "    'OP_2'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_1'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_2'  9  9   2   2 'OPEN' 1*   46.825   0.311  4332.346 1*  1*  'X'  22.123 / \n"
            " 'OP_1'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_1' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"
            "WCONINJE\n"
                "'OP_2' 'GAS' 'OPEN' 'RATE' 100 200 400 /\n"
            "/\n"

            "DATES             -- 2\n"
            " 20  JAN 2011 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_3'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_3'  9  9   1   1 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_3' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 3\n"
            " 15  JUN 2013 / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_2'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_1'  9  9   7  7 'SHUT' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"

            "DATES             -- 4\n"
            " 22  APR 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_4'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_4'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            " 'OP_3'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_4' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 5\n"
            " 30  AUG 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_5'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_5'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_5' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 6\n"
            " 15  SEP 2014 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_3' 'SHUT' 'ORAT' 20000  4* 1000 /\n"
            "/\n"

            "DATES             -- 7\n"
            " 9  OCT 2014 / \n"
            "/\n"
            "WELSPECS\n"
            "    'OP_6'       'OP'   9   9 1*     'OIL' 1*      1*  1*   1*  1*   1*  1*  / \n"
            "/\n"
            "COMPDAT\n"
            " 'OP_6'  9  9   3  9 'OPEN' 1*   32.948   0.311  3047.839 1*  1*  'X'  22.100 / \n"
            "/\n"
            "WCONPROD\n"
                "'OP_6' 'OPEN' 'ORAT' 20000  4* 1000 /\n"
            "/\n"
            "TSTEP            -- 8\n"
            "10 /"
            "/\n";

    const double i1BhpSet = 1.25;
    const double i2BhpSet = 2.35;
    const double i3BhpSet = 3.45;
    const double i1TemperatureSet = 4.55;
    const double i2TemperatureSet = 5.65;
    const double i3TemperatureSet = 6.75;
    const double i1WellRatesSet = 7.25;
    const double i2WellRatesSet = 8.35;
    const double i3WellRatesSet = 9.05;
    const double i1PerfRatesSet = 11.25;
    const double i2PerfRatesSet = 12.25;
    const double i3PerfRatesSet = 13.25;
    const double i1PerfPressSet = 14.25;
    const double i2PerfPressSet = 15.25;
    const double i3PerfPressSet = 16.25;

    DeckConstPtr deck = parser.parseString(input);
    EclipseStateConstPtr eclipseState = std::make_shared<const EclipseState>(deck);
    Opm::GridManager gridManager(deck);

    int numPhases = eclipseState->getNumPhases();
    size_t reportStep = 1;

    WellState wsResize;
    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 2, 9, numPhases);

    wsResize.bhp()[0] = i1BhpSet;
    wsResize.bhp()[1] = i2BhpSet;
    wsResize.temperature()[0] = i1TemperatureSet;
    wsResize.temperature()[1] = i2TemperatureSet;
    wsResize.wellRates()[0] = i1WellRatesSet;
    wsResize.wellRates()[5] = i2WellRatesSet;
    wsResize.perfRates()[0] =i1PerfRatesSet;
    wsResize.perfRates()[1] = i2PerfRatesSet;
    wsResize.perfRates()[7] = i3PerfRatesSet;
    wsResize.perfPress()[0] = i1PerfPressSet;
    wsResize.perfPress()[1] = i2PerfPressSet;
    wsResize.perfPress()[7] = i3PerfPressSet;

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 3, 10, numPhases);

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 3, 16, numPhases);

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 4, 30, numPhases);

    wsResize.bhp()[3] = i3BhpSet;
    wsResize.temperature()[3] = i3TemperatureSet;
    wsResize.wellRates()[9] = i3WellRatesSet;

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 5, 37, numPhases);

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 4, 29, numPhases);

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 5, 36, numPhases);

    wsResize.resize(*eclipseState, reportStep++);
    testSizes(wsResize, 5, 36, numPhases);

    double i1BhpGet = wsResize.bhp()[0];
    double i2BhpGet = wsResize.bhp()[1];
    double i3BhpGet = wsResize.bhp()[2];
    double i1TemperatureGet = wsResize.temperature()[0];
    double i2TemperatureGet = wsResize.temperature()[1];
    double i3TemperatureGet = wsResize.temperature()[2];
    double i1WellRatesGet = wsResize.wellRates()[0];
    double i2WellRatesGet = wsResize.wellRates()[5];
    double i3WellRatesGet = wsResize.wellRates()[6];
    double i1PerfRatesGet = wsResize.perfRates()[0];
    double i2PerfRatesGet = wsResize.perfRates()[1];
    double i3PerfRatesGet = wsResize.perfRates()[6];
    double i1PerfPressGet = wsResize.perfPress()[0];
    double i2PerfPressGet = wsResize.perfPress()[1];
    double i3PerfPressGet = wsResize.perfPress()[6];

    BOOST_CHECK_EQUAL(i1BhpGet, i1BhpSet);
    BOOST_CHECK_EQUAL(i2BhpGet, i2BhpSet);
    BOOST_CHECK_EQUAL(i3BhpGet, i3BhpSet);
    BOOST_CHECK_EQUAL(i1TemperatureGet, i1TemperatureSet);
    BOOST_CHECK_EQUAL(i2TemperatureGet, i2TemperatureSet);
    BOOST_CHECK_EQUAL(i3TemperatureGet, i3TemperatureSet);
    BOOST_CHECK_EQUAL(i1WellRatesGet, i1WellRatesSet);
    BOOST_CHECK_EQUAL(i2WellRatesGet, i2WellRatesSet);
    BOOST_CHECK_EQUAL(i3WellRatesGet, i3WellRatesSet);
    BOOST_CHECK_EQUAL(i1PerfRatesGet, i1PerfRatesSet);
    BOOST_CHECK_EQUAL(i2PerfRatesGet, i2PerfRatesSet);
    BOOST_CHECK_EQUAL(i3PerfRatesGet, i3PerfRatesSet);
    BOOST_CHECK_EQUAL(i1PerfPressGet, i1PerfPressSet);
    BOOST_CHECK_EQUAL(i2PerfPressGet, i2PerfPressSet);
    BOOST_CHECK_EQUAL(i3PerfPressGet, i3PerfPressSet);
}
