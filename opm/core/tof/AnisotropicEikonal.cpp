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

#include <opm/core/tof/AnisotropicEikonal.hpp>
#include <opm/core/grid/GridUtilities.hpp>
#include <opm/core/grid.h>
#include <set>

namespace Opm
{

    /// Construct solver.
    /// \param[in] grid      A 2d grid.
    AnisotropicEikonal2d::AnisotropicEikonal2d(const UnstructuredGrid& grid)
	: grid_(grid)
    {
	if (grid.dimensions != 2) {
	    OPM_THROW(std::logic_error, "Grid for AnisotropicEikonal2d must be 2d.");
	}
	cell_neighbours_ = vertexNeighbours(grid);
	orderCounterClockwise(grid, cell_neighbours_);
    }

    /// Solve the eikonal equation.
    /// \param[in]  metric            Array of metric tensors, M, for each cell.
    /// \param[in]  startcells        Array of cells where u = 0 at the centroid.
    /// \param[out] solution          Array of solution to the eikonal equation.
    void AnisotropicEikonal2d::solve(const double* metric,
				     const std::vector<int>& startcells,
				     std::vector<double>& solution)
    {
	// The algorithm used is described in J.A. Sethian and A. Vladimirsky,
	// "Ordered Upwind Methods for Static Hamilton-Jacobi Equations".
	// Notation in comments is as used in that paper: U is the solution,
	// and q is the boundary condition. One difference is that we talk about
	// grid cells instead of mesh points.
	//
	// Algorithm summary:
	// 1. Put all cells in Far. U_i = \inf.
	// 2. Move the startcells to Accepted. U_i = q(x_i)
	// 3. Move cells adjacent to startcells to Considered, evaluate
	//    U_i = min_{(x_j,x_k) \in NF(x_i)} G_{j,k}
	// 4. Find the Considered cell with the smallest value: r.
	// 5. Move cell r to Accepted. Update AcceptedFront.
	// 6. Move cells adjacent to r from Far to Considered.
	// 7. Recompute the value for all Considered cells within
	//    distance h * F_2/F1 from x_r. Use min of previous and new.
	// 8. If Considered is not empty, go to step 4.

	// 1. Put all cells in Far. U_i = \inf.
	const int num_cells = grid_.number_of_cells;
	const double inf = 1e100;
	solution.clear();
	solution.resize(num_cells, inf);
	is_accepted_.clear();
	is_accepted_.resize(num_cells, false);
        accepted_front_.clear();
	considered_.clear();
        considered_handles_.clear();
	is_considered_.clear();
	is_considered_.resize(num_cells, false);

	// 2. Move the startcells to Accepted. U_i = q(x_i)
	const int num_startcells = startcells.size();
	for (int ii = 0; ii < num_startcells; ++ii) {
	    is_accepted_[startcells[ii]] = true;
	    solution[startcells[ii]] = 0.0;
	}
	accepted_front_.insert(startcells.begin(), startcells.end());

	// 3. Move cells adjacent to startcells to Considered, evaluate
	//    U_i = min_{(x_j,x_k) \in NF(x_i)} G_{j,k}
	for (int ii = 0; ii < num_startcells; ++ii) {
	    const int scell = startcells[ii];
	    const int num_nb = cell_neighbours_[scell].size();
	    for (int nb = 0; nb < num_nb; ++nb) {
		const int nb_cell = cell_neighbours_[scell][nb];
		if (!is_accepted_[nb_cell] && !is_considered_[nb_cell]) {
		    const double value = computeValue(nb_cell, metric);
		    pushConsidered(std::make_pair(value, nb_cell));
		}
	    }
	}

        while (!considered_.empty()) {
            // 4. Find the Considered cell with the smallest value: r.
            const ValueAndCell r = topConsidered();

            // 5. Move cell r to Accepted. Update AcceptedFront.
            const int rcell = r.second;
            is_accepted_[rcell] = true;
            solution[rcell] = r.first;
            popConsidered();
            accepted_front_.insert(rcell);
            for (auto it = accepted_front_.begin(); it != accepted_front_.end();) {
                // Note that loop increment happens in the body of this loop.
                const int cell = *it;
                bool on_front = false;
                for (auto it2 = cell_neighbours_[cell].begin(); it2 != cell_neighbours_[cell].end(); ++it2) {
                    if (!is_accepted_[*it2]) {
                        on_front = true;
                        break;
                    }
                }
                if (!on_front) {
                    accepted_front_.erase(it++);
                } else {
                    ++it;
                }
            }

            // 6. Move cells adjacent to r from Far to Considered.
            for (auto it = cell_neighbours_[rcell].begin(); it != cell_neighbours_[rcell].end(); ++it) {
                const int nb_cell = *it;
                if (!is_accepted_[nb_cell] && !is_considered_[nb_cell]) {
                    assert(solution[nb_cell] == inf);
                    const double value = computeValue(nb_cell, metric);
                    pushConsidered(std::make_pair(value, nb_cell));
                }
            }

            // 7. Recompute the value for all Considered cells within
            //    distance h * F_2/F1 from x_r. Use min of previous and new.
            for (auto it = considered_.begin(); it != considered_.end(); ++it) {
                const int ccell = it->second;
                if (isClose(rcell, ccell, metric)) {
                    const double value = computeValue(ccell, metric);
                    if (value < it->first) {
                        // Update value for considered cell.
                        // Note that as solution values decrease, their
                        // goodness w.r.t. the heap comparator increase,
                        // therefore we may safely call the increase()
                        // modificator below.
                        considered_.increase(considered_handles_[ccell], std::make_pair(value, ccell));
                    }
                }
            }

            // 8. If Considered is not empty, go to step 4.
        }

    }





    bool AnisotropicEikonal2d::isClose(const int c1,
                                       const int c2,
                                       const double* metric) const
    {
        return true;
    }





    double AnisotropicEikonal2d::computeValue(const int cell,
                                              const double* metric) const
    {
	const auto& nbs = cell_neighbours_[cell];
	const int num_nbs = nbs.size();
        const double inf = 1e100;
	double val = inf;
	for (int ii = 0; ii < num_nbs; ++ii) {
	    const int n[2] = { nbs[ii], nbs[(ii+1) % num_nbs] };
            if (accepted_front_.count(n[0]) && accepted_front_.count(n[1])) {
                const double cand_val = computeFromTri(cell, n[0], n[1], metric);
                val = std::min(val, cand_val);
            }
	}
        if (val == inf) {
            // Failed to find two accepted front nodes adjacent to this,
            // so we go for a single-neighbour update.
            for (int ii = 0; ii < num_nbs; ++ii) {
                if (accepted_front_.count(nbs[ii])) {
                    const double cand_val = computeFromLine(cell, nbs[ii], metric);
                    val = std::min(val, cand_val);
                }
            }
        }
        assert(val != inf);
	return val;
    }





    double AnisotropicEikonal2d::computeFromLine(const int cell,
                                                 const int from,
                                                 const double* metric) const
    {
        assert(!is_accepted_[cell]);
        assert(is_accepted_[from]);
        return 2.0;
    }





    double AnisotropicEikonal2d::computeFromTri(const int cell,
                                                const int n0,
                                                const int n1,
                                                const double* metric) const
    {
        assert(!is_accepted_[cell]);
        assert(is_accepted_[n0]);
        assert(is_accepted_[n1]);
        return 1.0;
    }





    const AnisotropicEikonal2d::ValueAndCell& AnisotropicEikonal2d::topConsidered() const
    {
	return considered_.top();
    }





    void AnisotropicEikonal2d::pushConsidered(const ValueAndCell& vc)
    {
	HeapHandle h = considered_.push(vc);
        considered_handles_[vc.second] = h;
        is_considered_[vc.second] = true;
    }





    void AnisotropicEikonal2d::popConsidered()
    {
        is_considered_[considered_.top().second] = false;
        considered_handles_.erase(considered_.top().second);
	considered_.pop();
    }





} // namespace Opm
