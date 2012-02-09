/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TWOPHASETRANSPORT_HPP_INCLUDED
#define TWOPHASETRANSPORT_HPP_INCLUDED

namespace Opm
{
    class IncompPropertiesInterface;
}

struct UnstructuredGrid;

void twophasetransport(
    const double *porevolume,
    const double *source,
    double dt,
    struct UnstructuredGrid *grid,
    const Opm::IncompPropertiesInterface* props,
    const double *darcyflux,
    double *saturation);

namespace Opm
{

    void reorderTransportTwophase(const double *porevolume,
				  const double *source,
				  const double dt,
				  const UnstructuredGrid *grid,
				  const IncompPropertiesInterface* props,
				  const double *darcyflux,
				  double *saturation);
} // namespace Opm


#endif /* TWOPHASETRANSPORT_HPP_INCLUDED */