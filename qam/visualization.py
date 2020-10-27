import vtk

def process_data(input_file, color):
    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName(input_file)

    extractor = vtk.vtkMarchingCubes()
    extractor.SetInputConnection(reader.GetOutputPort())
    extractor.SetValue(0, 1)

    smooth_filter = vtk.vtkSmoothPolyDataFilter()
    smooth_filter.SetInputConnection(extractor.GetOutputPort());
    smooth_filter.SetNumberOfIterations(5);
    smooth_filter.SetRelaxationFactor(0.5);
    smooth_filter.FeatureEdgeSmoothingOff();
    smooth_filter.BoundarySmoothingOn();
    smooth_filter.Update()

    cleaned = vtk.vtkCleanPolyData()
    cleaned.SetInputData(smooth_filter.GetOutput())

    normal_generator = vtk.vtkPolyDataNormals()
    normal_generator.ComputePointNormalsOn();
    normal_generator.ComputeCellNormalsOn();
    normal_generator.SetInputConnection(cleaned.GetOutputPort());
    normal_generator.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normal_generator.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(color)

    return normal_generator, actor

def create_distance_map(tumor, ablation):
    distance_filter = vtk.vtkDistancePolyDataFilter()
    distance_filter.SignedDistanceOn()
    distance_filter.SetInputConnection(1, tumor.GetOutputPort() )
    distance_filter.SetInputConnection(0, ablation.GetOutputPort() )
    distance_filter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(distance_filter.GetOutputPort())
    mapper.SetScalarRange( distance_filter.GetOutput().GetPointData().GetScalars().GetRange()[0],
                           distance_filter.GetOutput().GetPointData().GetScalars().GetRange()[1])
    print('Min distance: ', distance_filter.GetOutput().GetPointData().GetScalars().GetRange()[0])
    print('Max distance: ', distance_filter.GetOutput().GetPointData().GetScalars().GetRange()[1])
    mapper.SetScalarRange(-5, 15)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # actor.GetProperty().SetOpacity(0.75)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(mapper.GetLookupTable())
    scalar_bar.SetTitle("Distance")
    scalar_bar.SetNumberOfLabels(3)

    hueLut = vtk.vtkLookupTable()
    hueLut.SetTableRange(0, 10)
    hueLut.SetHueRange(0, 0.4)
    hueLut.SetSaturationRange(1, 1)
    hueLut.SetValueRange(1, 1)
    hueLut.Build()

    mapper.SetLookupTable(hueLut)
    scalar_bar.SetLookupTable(hueLut)

    return actor, scalar_bar

def visualize_3d_margin(tumor_file, ablation_file, output_file):
    colors = vtk.vtkNamedColors()

    colors.SetColor("TumorColor", [255, 64, 64, 255])
    colors.SetColor("AblationColor", [64, 64, 255, 64])
    colors.SetColor("BkgColor", [0, 0, 0, 255])

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(render_window)

    tumor_mesh, tumor_actor = process_data(tumor_file, colors.GetColor3d("TumorColor"))
    ablation_mesh, ablation_actor = process_data(ablation_file, colors.GetColor3d("AblationColor"))
    distance_actor, scalar_bar = create_distance_map(tumor_mesh, ablation_mesh)

    # setup scene
    camera = vtk.vtkCamera()
    camera.SetViewUp(0, 0, -1)
    camera.SetPosition(0, -1, 0)
    camera.SetFocalPoint(0, 0, 0)
    camera.ComputeViewPlaneNormal()
    camera.Azimuth(30.0)
    camera.Elevation(30.0)

    renderer.AddActor(distance_actor)
    renderer.AddActor(scalar_bar)
    # renderer.AddActor(ablation_actor)
    # renderer.AddActor(tumor_actor)
    renderer.SetActiveCamera(camera)
    renderer.ResetCamera()
    camera.Dolly(1.5)

    renderer.SetBackground(colors.GetColor3d("BkgColor"))
    render_window.SetSize(640, 480)

    renderer.ResetCameraClippingRange()

    exporter = vtk.vtkVRMLExporter()
    exporter.SetFileName(str(output_file))
    exporter.SetRenderWindow(render_window)
    exporter.Write()

    iren.Initialize()
    iren.Start()
