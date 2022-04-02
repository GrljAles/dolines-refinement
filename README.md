# dolines-refinement
Karst dolines shape refinement.

The aim of the project is to refine the shapes of karst dolines (karst terrain depressions) with some advanced techniques descirbed below. 

PROBLEM:

The development of remote sensing methods, computers and spatial data processing software has also affected the field of geomorphology. All three allow for faster analysis of larger amounts of data, but not more accurate. The reason for this is the poor mathematical definition of some geomorphological forms. The project focuses on the delimitation of karst depressions. For this purpose, we have developed a new approach that analyzes its half-sections when determining the edge of a depression and also effectively detects the edges of depressions on slopes. The method allows obtaining more accurate morphographic and morphometric data on depressions. The method was developed and tested on five depressions of the Podgrajsko podolja (SW Slovenia), and the results were compared with the results of a method that uses hydrological modeling to delimit the depressions. We observe a significant improvement in the results obtained with the new method.


The results of the research were published in the Dela: https://revije.ff.uni-lj.si/Dela/article/view/9737 in Slovenian language, abstract and sumary are translated to english though.

Inputs for the script are:
    1. Polygon shapefile with polygons defining depressions rims delineated using simple hydrological filling of the input DEM, subtracting the original input DEM from it, reclassifying the result, and exporting to vector.

    2. Raster DEM used for delineating dolines at first stage.

Output is polygon shapefile with refined depression rims/edges that are correctly defined also for depresions on slopes.