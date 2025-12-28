import numpy as np # for graphene function
import multiprocessing as mp
import time
from multiprocessing import Process
from ase import Atoms
from ase.io import write,read
from ase import Atoms
from ase.build import bulk, surface, add_adsorbate
from ase.visualize import view
#import io, py3Dmol
from ase.io import write, read
#from IPython.display import HTML

#grphene
from ase.neighborlist import neighbor_list
import numpy as np
import os
from pathlib import Path

def zno_slab():
  #BUILDING 2D ZNO_SLAB

  zno = bulk('ZnO', 'wurtzite', a=3.25, c=5.20)
  zno_supercell = zno.repeat((3, 3, 3))
  atoms = zno_supercell

  pseudos = {
      'Zn': 'Zn.pbe-dnl-kjpaw_psl.1.0.0.UPF',
      'O' : 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
  }

  input_data = {
      'control': {
          'calculation': 'relax',
          'prefix': 'zno_slab',
          'outdir': './tmp',
          'pseudo_dir': './pseudo',
      },
      'system': {
          'ibrav': 0,
          'ecutwfc': 40.0,
          'ecutrho': 320.0,
          'occupations': 'smearing',  # bulk semiconductor
          'smearing': 'mp',
          'degauss': 0.02,
      },
      'electrons': {
          'conv_thr': 1.0e-5,
          'mixing_beta': 0.3,
          'electron_maxstep': 60,
      },
  }

  kpts = (9, 9, 4)   # 3D k-point mesh
  if not os.path.exists("input"):
        os.makedirs("input")

  file_path = "input/zno_slab_scf.in"

  write(file_path,
        atoms,
        format='espresso-in',
        pseudopotentials=pseudos,
        input_data=input_data,
        kpts=kpts)

  print("Wrote zno_bulk_scf.in")

def graphene():
  #build graphene sheet
  a_cc = 1.42   # C–C bond length (Å)

  cell = [
    [3*a_cc,            0.0,                 0.0],
    [0.0,      np.sqrt(3)*a_cc,              0.0],
    [0.0,              0.0,                 15.0]
  ]

  positions = [
    [0.0,                 0.0,                         0.0],
    [a_cc,                0.0,                         0.0],
    [1.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
    [2.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
  ]

  graphene_unit = Atoms('C4', positions=positions,
                      cell=cell, pbc=[True, True, False])

  # 6×6 supercell
  graphene = graphene_unit.repeat((6, 6, 1))

  pseudos = {
    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',   # same name as on the HPC
  }

  input_data = {
      'control': {
          'calculation': 'relax',
          'prefix': 'graphene',
          'outdir': './tmp',
          'pseudo_dir': './pseudo',        # we'll create this on HPC
      },
      'system': {
          'ecutwfc': 50.0,
          'ecutrho': 400.0,
          'occupations': 'smearing',
          'smearing': 'mp',
          'degauss': 0.02,
	  'assume_isolated':'2D'
      },
      'electrons': {
          'conv_thr': 1.0e-8,
      },
  }

  kpts = (4, 4, 1)
  if not os.path.exists("input"):
        os.makedirs("input")

  file_path = "input/g_scf.in"
  write(file_path, graphene,
        format='espresso-in',
        pseudopotentials=pseudos,
        input_data=input_data,
        kpts=kpts)
  print("Wrote graphene_scf.in")

def MB():

  #BUILD MB DYE
  # 1. Read MB geometry from XYZ
  atoms = read('mb_dye.xyz')   # file you uploaded from Avogadro

  print("Loaded MB with", len(atoms), "atoms")
  print("Elements:", set(a.symbol for a in atoms))

  # 2. Put MB into a 20×20×20 Å box (isolated molecule)
  cell = [20.0, 20.0, 20.0]
  atoms.set_cell([[cell[0], 0.0,      0.0],
                      [0.0,      cell[1], 0.0],
                      [0.0,      0.0,     cell[2]]],
                    scale_atoms=False)
  atoms.center()   # center the molecule in the box

  # 3. Define pseudopotentials (must match the files in ./pseudo on the HPC)
  pseudos = {
          'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
          'H': 'H.pbe-rrkjus.UPF',
          'N': 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
          'S': 'S.pbe-n-kjpaw_psl.1.0.0.UPF',
          'Cl': 'Cl.pbe-n-kjpaw_psl.1.0.0.UPF',  # uncomment if your MB has Cl-
      }

  # 4. QE input parameters for SCF-only MB
  input_data = {
          'control': {
              'calculation': 'scf',
              'prefix': 'mb',
              'outdir': './tmp',
              'pseudo_dir': './pseudo',
          },
          'system': {
              'ibrav': 0,          # free cell; use CELL_PARAMETERS from atoms
              'ecutwfc': 40.0,
              'ecutrho': 320.0,
              'occupations': 'fixed',   # molecule
          },
          'electrons': {
              'conv_thr': 1.0e-5,
              'mixing_beta': 0.3,
              'electron_maxstep': 60,
          },
      }

  # 5. Gamma-only k-points for isolated molecule
  kpts = (1, 1, 1)

  # 6. Write QE input file
  if not os.path.exists("input"):
        os.makedirs("input")

  file_path = "input/mb_scf.in"

  write(file_path,
            atoms,
            format='espresso-in',
            pseudopotentials=pseudos,
            input_data=input_data,
            kpts=kpts)

  print("Wrote mb_scf.in")

def g_zno():
    a_cc = 1.42   # C–C bond length (Å)

    cell = [
        [3*a_cc,            0.0,                 0.0],
        [0.0,      np.sqrt(3)*a_cc,              0.0],
        [0.0,              0.0,                 15.0]
    ]

    positions = [
        [0.0,                 0.0,                         0.0],
        [a_cc,                0.0,                         0.0],
        [1.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
        [2.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
    ]

    graphene_unit = Atoms('C4', positions=positions,
                          cell=cell, pbc=[True, True, False])

    # 6×6 supercell
    graphene = graphene_unit.repeat((6, 6, 1))

    zno_bulk = bulk('ZnO', 'wurtzite', a=3.25, c=5.20)
    zno_slab = surface(zno_bulk, (0, 0, 1), layers=8, vacuum=15.0)
    zno_slab.center(axis=2)
    zno_slab_big = zno_slab.repeat((3, 3, 1))

    # Copy graphene so we don't modify the original
    slab = graphene.copy()

    # ZnO molecule to adsorb
    ads = zno_slab_big.copy()

    # Height (distance between lowest atom in ZnO and average graphene z)
    height = 1.8 # Å

    # Put it roughly at the center of the graphene cell
    x = slab.cell[0,0] / 2
    y = slab.cell[1,1] / 2

    add_adsorbate(slab, ads, height=height, position=(x, y))

    hybrid_cluster = slab
    hybrid_2x2 = hybrid_cluster.repeat((1, 1, 1))  # 2×2 tiling for visualization
    hybrid_2x2.center(vacuum=10.0, axis=2)

    pseudos = {
        'C':  'C.pbe-n-kjpaw_psl.1.0.0.UPF',
        'O':  'O.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Zn': 'Zn.pbe-dnl-kjpaw_psl.1.0.0.UPF',
    }

    input_data = {
        'control': {
            'calculation': 'relax',
            'prefix': 'gzn_hybrid',
            'outdir': './tmp',
            'pseudo_dir': './pseudo',
        },
        'system': {
            'ibrav': 0,
            'ecutwfc': 40.0,
            'ecutrho': 320.0,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02,
        },
        'electrons': {
            'conv_thr': 1.0e-5,
            'mixing_beta': 0.3,
            'electron_maxstep': 80,
        },
        'ions': {
            'ion_dynamics': 'bfgs',
        },
    }

    kpts = (1, 1, 1)

    if not os.path.exists("input"):
        os.makedirs("input")

    file_path = "input/g_zno_scf.in"
    write(file_path,
          hybrid_2x2,
          format='espresso-in',
          pseudopotentials=pseudos,
          input_data=input_data,
          kpts=kpts)
    print("Wrote gzn_hybrid_relax.in")

def mb_gzno():

    a_cc = 1.42   # C–C bond length (Å)

    cell = [
        [3*a_cc,            0.0,                 0.0],
        [0.0,      np.sqrt(3)*a_cc,              0.0],
        [0.0,              0.0,                 15.0]
    ]

    positions = [
        [0.0,                 0.0,                         0.0],
        [a_cc,                0.0,                         0.0],
        [1.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
        [2.5*a_cc,  np.sqrt(3)/2*a_cc,                    0.0],
    ]

    graphene_unit = Atoms('C4', positions=positions,
                          cell=cell, pbc=[True, True, False])

    # 6×6 supercell
    graphene = graphene_unit.repeat((6, 6, 1))

    zno_bulk = bulk('ZnO', 'wurtzite', a=3.25, c=5.20)
    zno_slab = surface(zno_bulk, (0, 0, 1), layers=8, vacuum=15.0)
    zno_slab.center(axis=2)
    zno_slab_big = zno_slab.repeat((3, 3, 1))

    # Copy graphene so we don't modify the original
    slab = graphene.copy()

    # ZnO molecule to adsorb
    ads = zno_slab_big.copy()

    # Height (distance between lowest atom in ZnO and average graphene z)
    height = 1.8 # Å

    # Put it roughly at the center of the graphene cell
    x = slab.cell[0,0] / 2
    y = slab.cell[1,1] / 2

    add_adsorbate(slab, ads, height=height, position=(x, y))

    hybrid_cluster = slab
    hybrid_2x2 = hybrid_cluster.repeat((1, 1, 1))  # 2×2 tiling for visualization
    hybrid_2x2.center(vacuum=10.0, axis=2)

    ads = read('mb_dye.xyz')
    system = hybrid_2x2.copy()
    height = 3.5              # Å above graphene as a starting guess
    x = system.cell[0, 0] / 2
    y = system.cell[1, 1] / 2

    add_adsorbate(system, ads, height=height, position=(x, y))
    mb_on_gz = system
    pseudos = {
        'C':  'C.pbe-n-kjpaw_psl.1.0.0.UPF',
        'H':  'H.pbe-rrkjus.UPF',
        'N':  'N.pbe-n-kjpaw_psl.1.0.0.UPF',
        'S':  'S.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Cl': 'Cl.pbe-n-kjpaw_psl.1.0.0.UPF',   # if your MB has Cl-
        'O':  'O.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Zn': 'Zn.pbe-dnl-kjpaw_psl.1.0.0.UPF', # or Zn.pbe.UPF if you renamed
    }

    # 4. QE input parameters
    input_data = {
        'control': {
            'calculation': 'relax',     # relax MB+G/ZnO; use 'scf' if geometry is already relaxed
            'prefix': 'mb_gz',
            'outdir': './tmp',
            'pseudo_dir': './pseudo',
        },
        'system': {
            'ibrav': 0,
            'ecutwfc': 40.0,
            'ecutrho': 320.0,
            'occupations': 'smearing',  # safer for G/ZnO + MB interface
            'smearing': 'mp',
            'degauss': 0.02,
        },
        'electrons': {
            'conv_thr': 1.0e-5,
            'mixing_beta': 0.3,
            'electron_maxstep': 80,
        },
        'ions': {
            'ion_dynamics': 'bfgs',
        },
    }

    # 5. k-point mesh: small for first test (increase later)
    kpts = (1, 1, 1)

    if not os.path.exists("input"):
        os.makedirs("input")

    file_path = "input/mb_gzno.in"

    # 6. Write QE input
    write(file_path,
          mb_on_gz,
          format='espresso-in',
          pseudopotentials=pseudos,
          input_data=input_data,
          kpts=kpts)

    print("Wrote mb_gz_relax.in")



if __name__ == "__main__":
    # Create processes
    p1 = Process(target=zno_slab)
    p2 = Process(target=graphene)
    p3 = Process(target=MB)
    p4 = Process(target=g_zno)
    p5 = Process(target=mb_gzno)
    # Start processes (parallel execution)
    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()


    # Wait for completion
    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()

    print("All molecules built in parallel.")
