import anuga, time, os
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain
from anuga import Inlet_operator
from anuga import myid, distribute, finalize, barrier
from anuga.operators.rate_operators import Polygonal_rate_operator
import anuga.utilities.quantity_setting_functions as qs
import anuga.utilities.spatialInputUtil as su
import anuga.utilities.plot_utils as util
import json
import logging
from osgeo import ogr, gdal, osr

if __name__ != '__main__':
    from ..database_utils import write_percentage_complete


def start_sim(run_id, Runs, scenario_name, Scenario, session, **kwargs):
    yieldstep = kwargs['yieldstep']
    finaltime = kwargs['finaltime']
    logger = logging.getLogger(run_id)
    max_triangle_area = kwargs['max_triangle_area']
    logger.info('Starting hydrata_project')

    if run_id == 'local_run':
        base_dir = os.getcwd()
    else:
        base_dir = os.getcwd() + '/base_dir/%s/' % run_id

    outname = run_id
    meshname = base_dir + 'outputs/' + run_id + '.msh'

    def get_filename(data_type, file_type):
        files = os.listdir('%sinputs/%s' % (base_dir, data_type))
        filename = '%sinputs/%s/%s' % (base_dir,
                                       data_type,
                                       [f for f in files if f[-4:] == file_type][0]
                                       )
        return filename

    boundary_data_filename = get_filename('boundary_data', '.shp')
    elevation_data_filename = get_filename('elevation_data', '.tif')
    try:
        structures_filename = get_filename('structures', '.shp')
    except OSError as e:
        structures_filename = None
    try:
        rain_data_filename = get_filename('rain_data', '.shp')
    except OSError as e:
        rain_data_filename = None
    try:
        inflow_data_filename = get_filename('inflow_data', '.shp')
    except OSError as e:
        inflow_data_filename = None
    try:
        friction_data_filename = get_filename('friction_data', '.shp')
    except OSError as e:
        friction_data_filename = None

    logger.info('boundary_data_filename: %s' % boundary_data_filename)
    logger.info('structures_filename: %s' % structures_filename)
    logger.info('rain_data_filename: %s' % rain_data_filename)
    logger.info('inflow_data_filename: %s' % inflow_data_filename)
    logger.info('friction_data_filename: %s' % friction_data_filename)
    logger.info('elevation_data_filename: %s' % elevation_data_filename)

    # create a list of project files
    vector_filenames = [
        boundary_data_filename,
        structures_filename,
        rain_data_filename,
        inflow_data_filename,
        friction_data_filename
    ]

    # set the projection system for ANUGA calculations from the geotiff elevation data
    elevation_data_gdal = gdal.Open(elevation_data_filename)
    project_spatial_ref = osr.SpatialReference()
    project_spatial_ref.ImportFromWkt(elevation_data_gdal.GetProjectionRef())
    project_spatial_ref_epsg_code = int(project_spatial_ref.GetAttrValue("AUTHORITY", 1))

    # check the spatial reference system of the project files matches that of the calculation
    for filename in vector_filenames:
        if filename:
            prj_text = open(filename[:-4] + '.prj').read()
            srs = osr.SpatialReference()
            srs.ImportFromESRI([prj_text])
            srs.AutoIdentifyEPSG()
            logger.info('filename is: %s' % filename)
            logger.info('EPSG is: %s' % srs.GetAuthorityCode(None))
            if str(srs.GetAuthorityCode(None)) != str(project_spatial_ref_epsg_code):
                logger.warning('warning spatial refs are not maching: %s, %s' % (
                    srs.GetAuthorityCode(None),
                    project_spatial_ref_epsg_code
                ))

    logger.info('Setting up structures...')
    if structures_filename:
        structures = []
        logger.info('processing structures from :%s' % structures_filename)
        ogr_shapefile = ogr.Open(structures_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            structure = json.loads(ogr_layer_feature.GetGeometryRef().ExportToJson())['coordinates'][0]
            structures.append(structure)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

        logger.info('structures: %s' % structures)
    else:
        logger.warning('warning: no structures found.')
        structures = None

    logger.info('Setting up friction...')
    frictions = []
    if friction_data_filename:
        logger.info('processing frictions from :%s' % friction_data_filename)
        ogr_shapefile = ogr.Open(friction_data_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            friction_poly = json.loads(ogr_layer_feature.GetGeometryRef().ExportToJson())['coordinates'][0]
            friction_value = float(ogr_layer_feature.GetField('mannings'))
            friction_couple = [friction_poly, friction_value]
            frictions.append(friction_couple)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

        frictions.append(['All', 0.04])
        logger.info('frictions: %s' % frictions)
    else:
        frictions.append(['All', 0.04])
        logger.info('warning: no frictions found.')

    logger.info('Setting up boundary conditions...')
    ogr_shapefile = ogr.Open(boundary_data_filename)
    ogr_layer = ogr_shapefile.GetLayer(0)
    ogr_layer_definition = ogr_layer.GetLayerDefn()
    logger.info('ogr_layer_definition.GetGeomType: %s' % ogr_layer_definition.GetGeomType())
    boundary_tag_index = 0
    bdy_tags = {}
    bdy = {}

    ogr_layer_feature = ogr_layer.GetNextFeature()
    while ogr_layer_feature:
        boundary_tag_key = ogr_layer_feature.GetField('bdy_tag_k')
        boundary_tag_value = ogr_layer_feature.GetField('bdy_tag_v')
        bdy_tags[boundary_tag_key] = [boundary_tag_index * 2, boundary_tag_index * 2 + 1]
        bdy[boundary_tag_key] = boundary_tag_value
        geom = ogr_layer_feature.GetGeometryRef().GetPoints()
        ogr_layer_feature = None
        ogr_layer_feature = ogr_layer.GetNextFeature()
        boundary_tag_index = boundary_tag_index + 1
        logger.info('bdy_tags: %s' % bdy_tags)
    logger.info('bdy: %s' % bdy)

    boundary_data = su.read_polygon(boundary_data_filename)

    create_mesh_from_regions(
        boundary_data,
        boundary_tags=bdy_tags,
        maximum_triangle_area=max_triangle_area,
        interior_regions=None,
        interior_holes=structures,
        filename=meshname,
        use_cache=False,
        verbose=True
    )

    domain = Domain(meshname, use_cache=False, verbose=True)
    domain.set_name(outname)
    domain.set_datadir(base_dir + '/outputs')
    logger.info(domain.statistics())
    poly_fun_pairs = [['Extent', elevation_data_filename.encode("utf-8")]]
    topography_function = qs.composite_quantity_setting_function(
        poly_fun_pairs,
        domain,
        nan_treatment='exception',
    )
    friction_function = qs.composite_quantity_setting_function(
        frictions,
        domain
    )
    domain.set_quantity('friction', friction_function, verbose=True)
    domain.set_quantity('stage', 0.0)
    domain.set_quantity('elevation', topography_function, verbose=True, alpha=0.99)
    domain.set_minimum_storable_height(0.005)

    logger.info('Applying rainfall...')
    if rain_data_filename:
        ogr_shapefile = ogr.Open(rain_data_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        rainfall = 0
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            rainfall = float(ogr_layer_feature.GetField('rate_mm_hr'))
            polygon = su.read_polygon(rain_data_filename)
            logger.info("applying Polygonal_rate_operator with rate, polygon:")
            logger.info(rainfall)
            logger.info(polygon)
            Polygonal_rate_operator(domain, rate=rainfall, factor=1.0e-6, polygon=polygon, default_rate=0.0)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

    logger.info('Applying surface inflows...')
    if inflow_data_filename:
        ogr_shapefile = ogr.Open(inflow_data_filename)
        ogr_layer = ogr_shapefile.GetLayer(0)
        ogr_layer_definition = ogr_layer.GetLayerDefn()
        ogr_layer_feature = ogr_layer.GetNextFeature()
        while ogr_layer_feature:
            in_fixed = float(ogr_layer_feature.GetField('in_fixed'))
            line = ogr_layer_feature.GetGeometryRef().GetPoints()
            logger.info("applying Inlet_operator with line, in_fixed:")
            logger.info(line)
            logger.info(in_fixed)
            Inlet_operator(domain, line, in_fixed, verbose=False)
            ogr_layer_feature = None
            ogr_layer_feature = ogr_layer.GetNextFeature()

    logger.info('Applying Boundary Conditions...')
    logger.info('Available boundary tags: %s' % domain.get_boundary_tags())

    Br = anuga.Reflective_boundary(domain)
    Bd = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
    Bt = anuga.Transmissive_boundary(domain)

    for key, value in bdy.iteritems():
        if value == 'Br':
            bdy[key] = Br
        elif value == 'Bd':
            bdy[key] = Bd
        elif value == 'Bt':
            bdy[key] = Bt
        else:
            logger.info('No matching boundary condition exists - please check your shapefile attributes in: %s' % boundary_data_filename)

    # set a default value for exterior & interior boundary if it is not already set
    try:
        bdy['exterior']
    except KeyError:
        bdy['exterior'] = Br
    try:
        bdy['interior']
    except KeyError:
        bdy['interior'] = Br

    logger.info('bdy: %s' % bdy)

    domain.set_boundary(bdy)

    domain = distribute(domain)
    logger.info('Beginning evolve phase...')
    for t in domain.evolve(yieldstep, finaltime):
        domain.write_time()
        print domain.timestepping_statistics()
        logger.info(domain.timestepping_statistics(track_speeds=True))
        percentage_complete = round(domain.time/domain.finaltime, 3)*100
        logger.info('%s percent complete' % percentage_complete)
        if run_id != 'local_run':
            write_percentage_complete(run_id, Runs, scenario_name, Scenario, session, percentage_complete)
    domain.sww_merge(delete_old=True)
    barrier()
    finalize()
    sww_file = base_dir + '/outputs/' + run_id + '.sww'
    sww_file = sww_file.encode('utf-8', 'ignore')  # sometimes run_id gets turned to a unicode object by celery
    util.Make_Geotif(
        swwFile=sww_file,
        output_quantities=['depth', 'velocity'],
        myTimeStep='max',
        CellSize=max_triangle_area,
        lower_left=None,
        upper_right=None,
        EPSG_CODE=project_spatial_ref_epsg_code,
        proj4string=None,
        velocity_extrapolation=True,
        min_allowed_height=1.0e-05,
        output_dir=(base_dir + '/outputs/'),
        bounding_polygon=boundary_data,
        internal_holes=structures,
        verbose=False,
        k_nearest_neighbours=3,
        creation_options=[]
    )
    logger.info("Done. Nice work.")

if __name__ == "__main__":
    # TODO: parse argv for local development
    start_sim('local_run', Runs='Runs', session='local_session', Scenario='Scenario', scenario_name='local_scenario')
