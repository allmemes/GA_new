import arcpy
import json

def create_rasters_for_inbound_check(dataBase, inputPolygon, cell_size, bufferDistance, outputRasterName):
    arcpy.env.workspace = dataBase
    polygonName = inputPolygon.replace("\\", "/").split("/")[-1].split(".")[0]

    # Load the polygon into geodatabase if not exists
    walk = arcpy.da.Walk(dataBase, datatype="FeatureClass")
    for i in walk:
        if polygonName not in i[2]:
            arcpy.FeatureClassToGeodatabase_conversion([inputPolygon], dataBase)

    try:
        print("Generating required rasters...")
        # 1, Create the line shapefile for the boundary polygon first.
        arcpy.management.PolygonToLine(polygonName, "boundaryLine")

        # 2, Create the outward buffer to prevent points from mutating out of extent
        arcpy.analysis.Buffer("boundaryLine", "outwardBuffer", str(bufferDistance) + " Meters", "LEFT")
        # Add a field and make it 1 for conversion.
        arcpy.management.AddField("outwardBuffer", "rasterID", "SHORT", 0)
        arcpy.management.CalculateField("outwardBuffer", "rasterID", 1, "PYTHON3")
        # Convert it to raster.
        arcpy.conversion.PolygonToRaster("outwardBuffer", "rasterID", "bufferRaster", "CELL_CENTER", "NONE", cell_size, "BUILD")
        # Reclassiify
        buffer_raster = arcpy.sa.Reclassify("bufferRaster", "Value", "1 0;NODATA 0", "DATA")
        newExtent = buffer_raster.extent
        buffer_raster.save("buffer_raster")

        # 3, Create the polygon raster.
        arcpy.management.AddField(polygonName, "rasterID", "SHORT", 0)
        arcpy.management.CalculateField(polygonName, "rasterID", 1, "PYTHON3")
        with arcpy.EnvManager(snapRaster="buffer_raster", extent = newExtent):
            arcpy.conversion.PolygonToRaster(polygonName, "rasterID", "polygonRaster", "CELL_CENTER", "NONE", cell_size, "BUILD")
        # Reclassify
        polygon_raster = arcpy.sa.Reclassify("polygonRaster", "Value", "1 1;NODATA 0", "DATA")
        # polygon_raster.save("polygon_raster")

        # 4, Buffer inward to obtain the area right on boundary.
        arcpy.analysis.Buffer("boundaryLine", "inwardBuffer", str(cell_size) + " Meters", "RIGHT")
        arcpy.management.AddField("inwardBuffer", "rasterID", "SHORT", 0)
        arcpy.management.CalculateField("inwardBuffer", "rasterID", 1, "PYTHON3")
        # Convert it to raster
        with arcpy.EnvManager(snapRaster="buffer_raster", extent = newExtent):
            arcpy.conversion.PolygonToRaster("inwardBuffer", "rasterID", "boundaryRaster", "CELL_CENTER", "NONE", cell_size, "BUILD")
        # Reclassify
        boundary_raster = arcpy.sa.Reclassify("boundaryRaster", "Value", "1 1;NODATA 0", "DATA")
        # boundary_raster.save("boundary_raster")

        # 5, Add all rasters up
        result_raster = polygon_raster + boundary_raster + buffer_raster
        result_raster.save(outputRasterName)
        print("Final raster created")
        return result_raster
    except:
        print("Required rasters already exist")
        return arcpy.Raster(outputRasterName)

def get_boundary_index(inRas):
    boundaryList = []
    arr = arcpy.RasterToNumPyArray(inRas, nodata_to_value=0)
    for i, j in enumerate(arr):
        for k, l in enumerate(j):
            if l == 2:
                boundaryList.append((i, k))
    return boundaryList

def search_min_distance(boundaryList, row, col):
    smallestDistance = float("inf")
    smallestRow = 0
    smallestCol = 0
    for i in boundaryList:
        distance = (i[0] - row)**2 + (i[1] - col)**2 
        if distance < smallestDistance:
            smallestDistance = distance
            smallestRow = i[0]
            smallestCol = i[1]
    return (smallestRow, smallestCol)

def write_raster_to_json(inRas, bonudaryList, outputRasterName):
    arr = arcpy.RasterToNumPyArray(inRas, nodata_to_value=0)
    print("Writing file to json")
    json_dict = {"extent": {"llx": inRas.extent.XMin,
                            "lly": inRas.extent.YMin,
                            "urx": inRas.extent.XMax,
                            "ury": inRas.extent.YMax}}
    for i, j in enumerate(arr):
        newList = j.tolist()
        for k, l in enumerate(newList):
            if l == 0:
                indice = search_min_distance(bonudaryList, i, k)
                newList[k] = (inRas.extent.XMin + indice[1] * cell_size + 1/2 * cell_size, inRas.extent.YMax - indice[0] * cell_size - 1/2 * cell_size)
        json_dict[i] = newList
    with open(outputRasterName + ".json", "w") as file:
        json.dump(json_dict, file)
    print("Finished")

def delete_raster(dataBase, rasterName):
    arcpy.env.workspace = dataBase
    rasters = arcpy.ListRasters("*", "GRID")
    for i in rasters:
        if i == rasterName:
            arcpy.Delete_management(i)
            print("Raster deleted")
            return True
    print("No such file")
    return False

if __name__ == "__main__":
    # Demonstration: format of the input parameters. The arcgis related databases needs to be created in advance.
    dataBase = r"\\engin-labs.m.storage.umich.edu\hangchen\windat.v2\Desktop\test\test.gdb"
    inputPolygon = "boundary"
    rasterName = "result_Raster"
    cell_size = 1

    # Note this buffferDistance is to avoid points being mutated to locations outisde the raster extent. This value
    # is related to the mutation_scale in the genetic algorithm, where mutation scale refers to the standard deviation.
    # This value refers to the maximum value that can be drawn from a normal distribution with 0 as mean and
    # mutation_scale as the standard deviation.
    bufferDistance = 50

    # preprocessing
    final_raster = create_rasters_for_inbound_check(dataBase, inputPolygon, cell_size, bufferDistance, rasterName)
    boundList = get_boundary_index(final_raster)
    write_raster_to_json(final_raster, boundList, rasterName)