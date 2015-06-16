/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLSTATE_HEADER_INCLUDED
#define OPM_WELLSTATE_HEADER_INCLUDED

#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <vector>
#include <cassert>
#include <cstddef>
#include <utility>
#include <iostream>
#include <string>

namespace Opm
{
    ///
    /// \brief The WellIndex struct
    ///     wellNumber
    ///         The index of the well as it occurs in the Schedule.OpenWells object.
    ///     wellCompletionOffset
    ///         Absolute offset, i.e. from the start, for the completions for the well.In effect, it is the index where this
    ///         well's completions will start. Used for the perfrates_ and perfpress_ vectors.
    ///     completionMap
    ///         Offset relative to wellCompletionOffset, for this completion. Used for the perfrates_ and perfpress_ vectors.
    struct WellIndex
    {
        size_t wellNumber;
        size_t wellCompletionOffset;
        std::map<std::size_t, size_t> completionMap;
    };

    /// The state of a set of wells.
    class WellState
    {
    public:
        /// Allocate and initialize if wells is non-null.
        /// Also tries to give useful initial values to the bhp() and
        /// wellRates() fields, depending on controls.  The
        /// perfRates() field is filled with zero, and perfPress()
        /// with -1e100.
        template <class State>
        void init(const Wells* wells, const State& state)
        {
            if (wells) {
                const int nw = wells->number_of_wells;
                const int np = wells->number_of_phases;
                bhp_.resize(nw);
                thp_.resize(nw);
                temperature_.resize(nw, 273.15 + 20); // standard temperature for now
                wellrates_.resize(nw * np, 0.0);
                for (int w = 0; w < nw; ++w) {
                    assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                    const WellControls* ctrl = wells->ctrls[w];
                    if (well_controls_well_is_stopped(ctrl)) {
                        // Stopped well:
                        // 1. Rates: assign zero well rates.
                        for (int p = 0; p < np; ++p) {
                            wellrates_[np*w + p] = 0.0;
                        }
                        // 2. Bhp: assign bhp equal to bhp control, if
                        //    applicable, otherwise assign equal to
                        //    first perforation cell pressure.
                        if (well_controls_get_current_type(ctrl) == BHP) {
                            bhp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            const int first_cell = wells->well_cells[wells->well_connpos[w]];
                            bhp_[w] = state.pressure()[first_cell];
                        }
                        // 3. Thp: assign thp equal to thp control, if applicable,
                        //    otherwise assign equal to bhp value.
                        if (well_controls_get_current_type(ctrl) == THP) {
                            thp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            thp_[w] = bhp_[w];
                        }
                    } else {
                        // Open well:
                        // 1. Rates: initialize well rates to match controls
                        //    if type is SURFACE_RATE.  Otherwise, we
                        //    cannot set the correct value here, so we
                        //    assign a small rate with the correct
                        //    sign so that any logic depending on that
                        //    sign will work as expected.
                        if (well_controls_get_current_type(ctrl) == SURFACE_RATE) {
                            const double rate_target = well_controls_get_current_target(ctrl);
                            const double * distr = well_controls_get_current_distr( ctrl );
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = rate_target * distr[p];
                            }
                        } else {
                            const double small_rate = 1e-14;
                            const double sign = (wells->type[w] == INJECTOR) ? 1.0 : -1.0;
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = small_rate * sign;
                            }
                        }

                        // 2. Bhp: initialize bhp to be target pressure if
                        //    bhp-controlled well, otherwise set to a
                        //    little above or below (depending on if
                        //    the well is an injector or producer)
                        //    pressure in first perforation cell.
                                    if (well_controls_get_current_type(ctrl) == BHP) {
                                        bhp_[w] = well_controls_get_current_target( ctrl );
                                    } else {
                                        const int first_cell = wells->well_cells[wells->well_connpos[w]];
                                        const double safety_factor = (wells->type[w] == INJECTOR) ? 1.01 : 0.99;
                                        bhp_[w] = safety_factor*state.pressure()[first_cell];
                        }

                        // 3. Thp: assign thp equal to thp control, if applicable,
                        //    otherwise assign equal to bhp value.
                                    if (well_controls_get_current_type(ctrl) == THP) {
                                        thp_[w] = well_controls_get_current_target( ctrl );
                                    } else {
                                        thp_[w] = bhp_[w];
                        }
                        }
                }

                // The perforation rates and perforation pressures are
                // not expected to be consistent with bhp_ and wellrates_
                // after init().
                perfrates_.resize(wells->well_connpos[nw], 0.0);
                perfpress_.resize(wells->well_connpos[nw], -1e100);
            }
        }

        ///
        /// \brief resize
        ///     This method resizes the vectors bhp_, temperature_, wellrates_, perfrates_ and perfpress,
        ///     to the sizes requiered for the report step that is about to be simulated.The values that
        ///     are already in these vectors, as a result of the previous report step that was simulated,
        ///     are repositioned to their new indices - taking into account that both wells and completions
        ///     can be both opened and closed.
        /// \param state
        ///     The EclipseState object.
        /// \param reportStep
        ///     The number of the report step that is about to be simulated.
        void resize(const EclipseState& state, const size_t reportStep)
        {
            std::vector<WellPtr> wells = state.getSchedule()->getOpenWells(reportStep);
            const size_t numPhases = state.getNumPhases();
            const int numCompl = getNumCompletions(wells, reportStep);

            const size_t numWells = wells.size();

            const std::vector<double> bhpOld = bhp_;
            const std::vector<double> temperatureOld = temperature_;
            const std::vector<double> wellRatesOld = wellrates_;
            const std::vector<double> perfratesOld = perfrates_;
            const std::vector<double> perfpressOld = perfpress_;

            std::map<std::string, WellIndex> newMap = createWellIndexMap(state, wells, reportStep);

            const double default_temp = 273.15 + 20.0;
            const double default_perfpress = -1e100;
            const double default_rates = 0.0;

            bhp_.resize(numWells);
            temperature_.resize(numWells, default_temp);
            wellrates_.resize(numWells * numPhases, default_rates);
            perfrates_.resize(numCompl, default_rates);
            perfpress_.resize(numCompl, default_perfpress);

            for (auto oldIter = wellIndexMap_.begin(); oldIter != wellIndexMap_.end(); ++oldIter) {
                const std::string& wellName = (*oldIter).first;
                auto newIter = newMap.find(wellName);

                if (newIter != newMap.end()) {
                    const WellIndex& oldWi = (*oldIter).second;
                    const size_t oldIndex = oldWi.wellNumber;

                    const WellIndex& newWi = (*newIter).second;
                    const size_t newIndex = newWi.wellNumber;

                    bhp_[newIndex] = bhpOld[oldIndex];
                    temperature_[newIndex] = temperatureOld[oldIndex];

                    for (int i = 0; i < numPhases; ++i) {
                        wellrates_[newIndex * numPhases + i] = wellRatesOld[oldIndex * numPhases + i];
                    }

                    for (auto oldComplIter = oldWi.completionMap.begin(); oldComplIter != oldWi.completionMap.end(); ++oldComplIter) {
                        size_t completionId = (*oldComplIter).first;
                        auto newComplIter = newWi.completionMap.find(completionId);

                        if (newComplIter != newWi.completionMap.end()) {
                            const size_t oldComplIndex = (*oldComplIter).second;
                            const size_t newComplIndex = (*newComplIter).second;

                            const size_t oldWellCompletionOffset = oldWi.wellCompletionOffset;
                            const size_t newWellCompletionOffset = newWi.wellCompletionOffset;

                            perfrates_[newWellCompletionOffset + newComplIndex] = perfratesOld[oldWellCompletionOffset + oldComplIndex];
                            perfpress_[newWellCompletionOffset + newComplIndex] = perfpressOld[oldWellCompletionOffset + oldComplIndex];
                        }
                    }
                }
            }

            wellIndexMap_ = newMap;
        }

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        /// One thp pressure per well.
        std::vector<double>& thp() { return thp_; }
        const std::vector<double>& thp() const { return thp_; }

        /// One temperature per well.
        std::vector<double>& temperature() { return temperature_; }
        const std::vector<double>& temperature() const { return temperature_; }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        /// One rate per well connection.
        std::vector<double>& perfRates() { return perfrates_; }
        const std::vector<double>& perfRates() const { return perfrates_; }

        /// One pressure per well connection.
        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

        size_t getRestartBhpOffset() const {
            return 0;
        }

        size_t getRestartPerfPressOffset() const {
            return bhp_.size();
        }

        size_t getRestartPerfRatesOffset() const {
            return getRestartPerfPressOffset() + perfpress_.size();
        }

        size_t getRestartTemperatureOffset() const {
            return getRestartPerfRatesOffset() + perfrates_.size();
        }

        size_t getRestartWellRatesOffset() const {
            return getRestartTemperatureOffset() + temperature_.size();
        }

    private:
        std::vector<double> bhp_;
        std::vector<double> thp_;
        std::vector<double> temperature_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;

        ///
        /// \brief wellIndexMap_
        ///     This map is used in the resize method. It holds information between report steps.
        std::map<std::string, WellIndex> wellIndexMap_;

        int getWellMapIndex(const std::string wellName) const
        {
            int retVal = -1;
            auto iter = wellIndexMap_.find(wellName);

            if (iter != wellIndexMap_.end()) {
                const WellIndex& wi = (*iter).second;
                retVal = wi.wellNumber;
            }

            return retVal;
        }

        int getNumCompletions(const std::vector<WellPtr>& wells, const int reportStep) const
        {
            int numCompl = 0;

            for (size_t i = 0; i < wells.size(); ++i) {
                int numComplForWell = getNumCompletions(wells[i], reportStep);
                numCompl += numComplForWell;
            }

            return numCompl;
        }

        int getNumCompletions(const WellPtr& well, const int reportStep) const
        {
            int numComplForWell = 0;
            CompletionSetConstPtr complSet = well->getCompletions(reportStep);

            for (size_t i = 0; i < complSet->size(); ++i) {
                WellCompletion::StateEnum complState = complSet->get(i)->getState();

                if (complState == WellCompletion::StateEnum::OPEN) {
                    numComplForWell++;
                }
            }

            return numComplForWell;
        }   

        std::map<std::string, WellIndex> createWellIndexMap(const EclipseState& state, const std::vector<WellPtr>& wells, const int reportStep)
        {
            EclipseGridConstPtr grid = state.getEclipseGrid();
            size_t numWells = wells.size();
            std::map<std::string, WellIndex> newMap;

            size_t totalWellCompletionOffset = 0;
            for (size_t i = 0; i < numWells; ++i) {
                WellIndex wi;
                wi.wellNumber = i;
                wi.wellCompletionOffset = totalWellCompletionOffset;

                totalWellCompletionOffset += getNumCompletions(wells[i], reportStep);
                fillCompletionsMap(grid, wells[i], reportStep, wi);
                std::pair <std::string, WellIndex> p(wells[i]->name(), wi);
                newMap.insert(p);
            }

            return newMap;
        }

        void fillCompletionsMap(const EclipseGridConstPtr& grid, const WellPtr& well, const int reportStep, WellIndex& wi) const
        {
            CompletionSetConstPtr complSet = well->getCompletions(reportStep);
            size_t open_completions_count = 0;

            for (size_t idx = 0; idx < complSet->size(); ++idx) {
                CompletionConstPtr completion = complSet->get(idx);

                if (completion->getState() == WellCompletion::StateEnum::OPEN) {
                    int i = completion->getI();
                    int j = completion->getJ();
                    int k = completion->getK();

                    size_t globalIndex = grid->getGlobalIndex(i, j, k);
                    std::pair <size_t, size_t> p(globalIndex, open_completions_count);
                    wi.completionMap.insert(p);
                    ++open_completions_count;
                }
            }
        }
    };

} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
