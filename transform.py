import meshio
import sys

mesh_filename =  "test.inp"

mesh = meshio.read(mesh_filename)

output_vtk = "test.vtk"
mesh.write(output_vtk)