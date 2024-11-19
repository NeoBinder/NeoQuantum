

def openmm_model2geom(model, charge):
    from openmm import unit
    topology = model.topology
    positions = model.positions.value_in_unit(unit.angstroms)
    is_open_shelled = False
    geom = ""
    line = ' {:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n'
    total_elec = 0
    for i in range(topology.n_atoms):
        #  x, y, z = positions[i][0], positions[i][1], positions[i][2]
        element = topology.atom(i).element
        symbol = element.symbol
        total_elec += element.atomic_number
        geom += line.format(symbol, *positions[i])

    if total_elec % 2 != 0:
        total_elec += charge
        if total_elec % 2 != 0:
            is_open_shelled = True

    return geom, is_open_shelled
