import numpy as np
import gdal
from osgeo.gdalconst import GA_ReadOnly
from osgeo import ogr
import os
from math import cos, sin, radians, degrees, atan2
import scipy.ndimage
from scipy.signal import find_peaks
from sklearn.linear_model import LinearRegression


def openRasterOneBand(inRasterP):
    """
    Opens single band raster and returns its data and metadata.
    Paramerters:
        inRasterP - path to raster file.
     Returns:
        inArray - numpy array of input band.
        inCols - number of input data columns.
        inRows - number of input data rows.
        inBands - number of input data bands.
        inNoData - no data value of input raster.
        inExtension - extension of input raster.
        inDriver - gdal driver of input raster
        inGeotransform - list of six affine transform C describing the relationship between raster positions (in pixel/line coordinates) and georeferenced coordinates.
        inProj - georeferencing coordinate system of input data.
        inProjRef - projection coordinate system of the image in OpenGIS WKT format.
    """
    inRaster = gdal.Open(inRasterP, GA_ReadOnly)
    inBand = inRaster.GetRasterBand(1)
    inNoData = inBand.GetNoDataValue()
    inArray = inBand.ReadAsArray().astype("float32")
    inExtension = os.path.splitext(inRasterP)[1]
    inCols = inRaster.RasterXSize
    inRows = inRaster.RasterYSize
    inBands = inRaster.RasterCount
    inDriver = inRaster.GetDriver().ShortName
    inGeotransform = inRaster.GetGeoTransform()
    inProj = inRaster.GetProjection()
    inProjRef = inRaster.GetProjectionRef()
    del inRaster
    del inBand

    return (inArray, inCols, inRows, inBands, inNoData, inExtension, inDriver, inGeotransform, inProj, inProjRef)


def moving_average(a, n=3):
    """
    Calculates a moving average over 1d array and returns smoothed array.
    :param a: 1d array
    :param n: Size of the neighbourhood to consider while averaging.
    :return: Smoothed 1d array.
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]

    return ret[n - 1:] / n


def point_on_circle(center, angle, radius):
    """
    Finding the x,y coordinates on circle, based on given angle and radius of the circle.
    :param center: Coordinates of the center of the circle.
    :param angle: Angle in degrees
    :param radius: Distance from center or radius of the circle.
    :return: X and Y coordinate of the point on the circle.
    """
    angle = radians(angle)

    x = center[0] + (radius * cos(angle))
    y = center[1] + (radius * sin(angle))

    return x, y


def extendedProfile(p0, p1, extension, surface):
    """
    Calculates a straight profile between two points over a given surface.
    :param p0: Coordinate pair of first point defining the profile, given as a list of y and x
    :param p1: Coordinate pair of second point defining the profile, given as a list of  y and x
    :param extension: Factor for profile length multiplication. Given as a positive number.
    :param surface: 2d numpy array to calculate profile over.
    :return: 1d array of surface values along profile, porfile length and it"s angle.
    """
    # Ger single coordinates from pairs.
    y0 = p0[0]
    x0 = p0[1]
    y1 = p1[0]
    x1 = p1[1]

    # Distance between p0 and p1.
    length = int(np.hypot(y1 - y0, x1 - x0))

    # Get azimuth between p0 and p1
    angle = degrees(atan2((x1 - x0), (y1 - y0)))

    # Calculate new point coordinates using extnesion distance in azimuth direction.
    y2, x2 = point_on_circle([y0, x0], angle, length * extension)

    # Coordinate pairs along profile to extract DEM values at.
    x, y = np.linspace(
        x0, x2, length * extension), np.linspace(y0, y2, length * extension)

    # Extract the values along the line
    profileValues = scipy.ndimage.map_coordinates(surface, np.vstack((x, y)))

    return profileValues, length, angle


def linearProjection(k, kx, n):
    """
    Calculates an estimation of y values for x using the linear regression parameters.
    :param k: Slope coefficient of the regression line
    :param kx: x value
    :param n: y intercept of the regression line.
    :return: A y value for given x.
    """
    return k * kx + n


def regressionLine(xData, yData):
    """
    Calculates a linear regression parameters and returns a values array of the line of best fit estimate.
    :param xData: Independant variable as numpy array.
    :param yData: Dependant variable as numpy array.
    :return: A numpy array with the linear projection of x and y data.
    """
    # Fit x and y data
    linModel = LinearRegression().fit(xData, yData)
    # Calcualte linear projection of every pair.
    regLin = linearProjection(linModel.coef_[0], xData, linModel.intercept_)
    # reshape to get a 1d array.
    regLin = regLin.reshape(1, -1)[0]

    return regLin


# Create results shapefile.
driver = ogr.GetDriverByName("Esri Shapefile")
# <= Define full path to putput shapefile here.
resDataset = driver.CreateDataSource("path/to/result/shapefile.shp")
resLayer = resDataset.CreateLayer("", None, ogr.wkbPolygon)
# Create the id field.
resLayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

# Define inputs.
# <= Define full path to input DEM that will be used for evelation inforamtion extraction.
raster = openRasterOneBand("path/to/input/DEM.tif")
# <= Define full path to input polygon shapefile defining the basic shapes of the terrain depressions.
depressionDataset = ogr.Open("path/to/input/terrain/depressions/shapefile.shp")

# Get raster as numpy array and iterate depressions features.
z = raster[0]

# Fix noData values to raster average.
z = np.where(z == raster[4], np.mean(z[z > raster[4]]), z)

# Get the upper left cornner coordinate of the input raster.
ulx, uly = raster[7][3], raster[7][0]

layer = depressionDataset.GetLayer()

id = 0
for feature in layer:
    geom = feature.GetGeometryRef()

    # Calculate the minimum bounding geomerty (convex hull) around polygon to simplify its geometry and remove possible comlications caused by its jagged rim.
    geom.ConvexHull()
    # Provide polygon rim with sufficient vertices to calculate as many profiles as possible.
    geom.Segmentize(1)
    geom.CloseRings()

    # Get depression center point coordinates and convert it to the image pixel coordinates using raster upper left corner coorddinates.
    cy, cx, cz = geom.Centroid().GetPoint()
    x0, y0 = np.abs(cx - ulx), np.abs(cy - uly)

    # Geometry for storing the results
    peakRing = ogr.Geometry(ogr.wkbLinearRing)

    # Iterate over each vertex of depressions rim.
    for depressioni in geom:
        for p in range(0, depressioni.GetPointCount()):
            # Get the coordinates and convert to pixel coordinates.
            y1, x1, z1 = depressioni.GetPoint(p)
            y1, x1 = np.abs(y1 - uly), np.abs(x1 - ulx)

            # Extract the profile values.
            zi, length, angle = extendedProfile([y0, x0], [y1, x1], 4, z)

            # Smoothen the profile and calculate its slope (first derivative).
            zi = moving_average(zi, n=3)
            mi = np.diff(zi) / np.diff(range(len(zi)))

            # Distribute length values along profile.
            xi = np.arange(len(zi)).reshape((-1, 1))
            # Calculate the linear regression parameters over the profile.
            rl = regressionLine(xi, zi)

            """
            Calculate the difference between the actual values over the profile and its linear regression to sort of
            normalize the profile giving the areas below the line negative values and above positive. This transforms 
            the profile so the areas above the line appear as peaks. Its maximums can be detected as rim of the actual depression.
            """
            mi = zi - rl

            """
            Find the peaks in then profile. First search if there is a peak in the elevation profile. If no peak can be 
            determined from it use the difference between the profile and its linear regression. If this fails too,
            start reducing the prominence parameter by dividing it by two for every iteration that does not identifies 
            peak. If also this fails the depression rim has no solution and the profile in question is skipped.  
            """
            prom = 10
            peaks2, _ = find_peaks(mi, prominence=prom)
            if len(peaks2) > 0:
                py, px = point_on_circle([y0, x0], angle, peaks2[0])
                peakRing.AddPoint(uly + py, ulx - px)
            else:
                peaks2, _ = find_peaks(mi, prominence=prom)
                if len(peaks2) > 0:
                    py, px = point_on_circle([y0, x0], angle, peaks2[0])
                    peakRing.AddPoint(uly + py, ulx - px)
                while len(peaks2) == 0:
                    prom = prom / 2
                    peaks2, _ = find_peaks(mi, prominence=prom)
                    if len(peaks2) > 0:
                        py, px = point_on_circle([y0, x0], angle, peaks2[0])
                        peakRing.AddPoint(uly + py, ulx - px)
                    if prom == 0.:
                        break

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(peakRing)
    poly = poly.ConvexHull()
    peakRing = None
    # Create the feature and set values
    defn = resLayer.GetLayerDefn()
    feat = ogr.Feature(defn)
    defn = None
    feat.SetField("id", id)
    feat.SetGeometry(poly)
    poly = None
    resLayer.CreateFeature(feat)
    feat = None
    id += 1
# close the shapefile.
resDataset.Destroy()
