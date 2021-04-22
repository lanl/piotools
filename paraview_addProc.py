input = self.GetInputDataObject(0, 0)
output = self.GetOutputDataObject(0)
output.CopyStructure(input)
array = vtk.vtkIntArray()
array.SetName("ProcID")
array.SetNumberOfComponents(1)
array.SetNumberOfTuples(input.GetNumberOfCells())
procID = numpy.load("/tmp/proc.npz")["procID"]
print(len(procID))
for i in range(input.GetNumberOfCells()):
    array.SetTuple1(i, procID[i])
output.GetCellData().AddArray(array)
print(input.GetNumberOfCells())
