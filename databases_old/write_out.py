from ase.db import connect


with connect('metals_dissolution.db') as db, connect('../../databases/metals_dissolution_out.db') as db2, connect('metals_dissolution_out.db') as db_out:
    for row in db.select():

        try:
            atoms = db2.get_atoms(metal=row.metal,surface=row.surface,defect=row.defect)
            db_out.write(atoms,idx=row.id,metal=row.metal,surface=row.surface,defect=row.defect)
        except KeyError:
            db_out.reserve(idx=row.id,metal=row.metal,surface=row.surface,defect=row.defect)