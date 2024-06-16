# -*- coding: utf-8 -*-
import copy
import math
TOL_ERROR = 0.0000001   # Error
from shapely.geometry import Polygon


def almost_equal(a, b, tolerance=None):
    """
    returns true if two points are approximately equal
    :param a: value
    :param b: value
    :param tolerance: Error value
    :return:
    """
    if tolerance is None:
        tolerance = TOL_ERROR
    return abs(a - b) < tolerance


def normalize_vector(v):
    """
    normalize vector into a unit vector
    :return:
    """
    if almost_equal(v['x'] * v['x'] + v['y'] * v['y'], 1):
        # given vector was already a unit vector
        return v
    inverse = 1
    if(math.sqrt(v['x']**2 + v['y']**2)!=0):
        inverse = 1.0 / math.sqrt(v['x']**2 + v['y']**2)

    return {'x': v['x']*inverse, 'y': v['y']*inverse}


def on_segment(A, B, p):
    """
    returns true if p lies on the line segment defined by AB, but not at any endpoints
    :param A:
    :param B:
    :param p:
    :return:
    """
    # vertical line
    if almost_equal(A['x'], B['x']) and almost_equal(p['x'], A['x']):
        if not almost_equal(p['y'], B['y']) and not almost_equal(p['y'], A['y']) and \
                        max(B['y'], A['y']) > p['y'] and p['y'] > min(B['y'], A['y']):
            return True
        else:
            return False
    # vertical line
    if almost_equal(A['y'], B['y']) and almost_equal(p['y'], A['y']):
        if not almost_equal(p['x'], B['x']) and not almost_equal(p['x'], A['x']) and \
                        max(B['x'], A['x']) > p['x'] and p['x'] > min(B['x'], A['x']):
            return True
        else:
            return False
    # range check
    if (p['x'] < A['x'] and p['x'] < B['x']) or (p['x'] > A['x'] and p['x'] > B['x']) or (
                    p['y'] < A['y'] and p['y'] < B['y']) or (p['y'] > A['y'] and p['y'] > B['y']):
        return False

    # exclude end points
    if (almost_equal(p['x'], A['x']) and almost_equal(p['y'], A['y'])) or (
                almost_equal(p['x'], B['x']) and almost_equal(p['y'], B['y'])):
        return False

    cross = (p['y'] - A['y']) * (B['x'] - A['x']) - (p['x'] - A['x']) * (B['y'] - A['y'])
    if abs(cross) > TOL_ERROR:
        return False
    dot = (p['x'] - A['x']) * (B['x'] - A['x']) + (p['y'] - A['y']) * (B['y'] - A['y'])
    if dot < 0 or almost_equal(dot, 0):
        return False

    len2 = (B['x'] - A['x']) * (B['x'] - A['x']) + (B['y'] - A['y']) * (B['y'] - A['y'])
    if dot > len2 or almost_equal(dot, len2):
        return False
    return True

def find_feasible_translation_vectors(A, B, touching):
    """
    generate translation vectors from touching vertices/edges
    returns feasible translation vectors
    """

    len_a = len(A['points'])
    len_b = len(B['points'])
    vectors = []
    for i in range(0, len(touching)):
        vertex_a = {'A': A['points'][touching[i]['A']], 'marked': True}

        prev_a_index = touching[i]['A'] - 1 
        prev_a_index = len_a - 1 if prev_a_index < 0 else prev_a_index  
        prev_a = A['points'][prev_a_index] 

        # adjacent B vertices
        vertex_b = {'A': B['points'][touching[i]['B']]} 
        prev_b_index = touching[i]['B'] - 1 
        next_b_index = touching[i]['B'] + 1 
        prev_b_index = len_b - 1 if prev_b_index < 0 else prev_b_index  
        next_b_index = 0 if next_b_index >= len_b else next_b_index  

        prev_b = B['points'][prev_b_index] 
        next_b = B['points'][next_b_index] 

        if touching[i]['type'] == 0:
            v_a1 = {
                'x': prev_a['x'] - vertex_a['A']['x'], 
                'y': prev_a['y'] - vertex_a['A']['y'], 
                'start': vertex_a['A'], 
                'end': prev_a  
            }

            v_b1 = {
                'x': vertex_b['A']['x'] - prev_b['x'], 
                'y': vertex_b['A']['y'] - prev_b['y'], 
                'start': vertex_b['A'], 
                'end': prev_b 
            }
            v_bb = {'start': {'x' : v_b1['start']['x'] + B['offsetx'], 'y' : v_b1['start']['y'] + B['offsety']}, 'end': { 'x' : v_b1['end']['x'] + B['offsetx'], 'y': v_b1['end']['y'] + B['offsety']}}
            num_vector = choose_translation_vector(v_a1, v_bb)

            v_a1, vector_intersaction_a = polygons_intersect_without_edge_touching(A, B, v_a1)
            v_b1, vector_intersaction_b = polygons_intersect_without_edge_touching(A, B, v_b1)

            if num_vector == 1:
                vectors.append(v_b1) if not vector_intersaction_b else None
            elif num_vector == 0:
                vectors.append(v_a1) if not vector_intersaction_a else None
            elif num_vector == 2:
                vectors.extend([v for v, intersects in [(v_b1, vector_intersaction_b), (v_a1, vector_intersaction_a)] if not intersects])
           

            v_b2 = {
                'x': vertex_b['A']['x'] - next_b['x'], 
                'y': vertex_b['A']['y'] - next_b['y'], 
                'start': next_b, 
                'end': vertex_b['A'] 
            }

            v_bb = {'start': {'x' : v_b2['start']['x'] + B['offsetx'], 'y' : v_b2['start']['y'] + B['offsety']}, 'end': { 'x' : v_b2['end']['x'] + B['offsetx'], 'y': v_b2['end']['y'] + B['offsety']}}
            num_vector = choose_translation_vector(v_a1, v_bb)

            v_b2, vector_intersaction_b = polygons_intersect_without_edge_touching(A, B, v_b2)

            if num_vector == 1:
                vectors.append(v_b2) if not vector_intersaction_b else None
            elif num_vector == 0:
                vectors.append(v_a1) if not vector_intersaction_a else None
            elif num_vector == 2:
                vectors.extend([v for v, intersects in [(v_b2, vector_intersaction_b), (v_a1, vector_intersaction_a)] if not intersects])


        
        elif touching[i]['type'] == 1:

            v = {
                'x': vertex_a['A']['x'] - (vertex_b['A']['x'] + B['offsetx']), 
                'y': vertex_a['A']['y'] - (vertex_b['A']['y'] + B['offsety']), 
                'start': prev_a, 
                'end': vertex_a['A'] 
            }
            v, vector_intersaction = polygons_intersect_without_edge_touching(A, B, v)
            if(not vector_intersaction):
                vectors.append(v)


            v = {
                'x': prev_a['x'] - (vertex_b['A']['x'] + B['offsetx']), 
                'y': prev_a['y'] - (vertex_b['A']['y'] + B['offsety']), 
                'start': vertex_a['A'], 
                'end': prev_a 
            }
            v, vector_intersaction = polygons_intersect_without_edge_touching(A, B, v)
            if(not vector_intersaction):
                vectors.append(v)

        elif touching[i]['type'] == 2:
            v = {
                'x': vertex_a['A']['x'] - (vertex_b['A']['x'] + B['offsetx']),
                'y': vertex_a['A']['y'] - (vertex_b['A']['y'] + B['offsety']),
                'start': prev_b,
                'end': vertex_b['A']
            }
            v, vector_intersaction = polygons_intersect_without_edge_touching(A, B, v)
            if(not vector_intersaction):
                vectors.append(v)

            v = {
                'x': vertex_a['A']['x'] - (prev_b['x'] + B['offsetx']),
                'y': vertex_a['A']['y'] - (prev_b['y'] + B['offsety']),
                'start': vertex_b['A'],
                'end': prev_b
            }
            v, vector_intersaction = polygons_intersect_without_edge_touching(A, B, v)
            if(not vector_intersaction):
                vectors.append(v)

    return vectors



def get_single_translation_vector(vectors, prevvector):
    """
    returns a single translation vector from the list of vectors
    """
    translate = None
    max_d = None

    for i in range(0, len(vectors)):
        if vectors[i]['x'] == 0 and vectors[i]['y'] == 0:
            continue
        
        # if this vector points us back to where we came from, ignore it.
        # ie cross product = 0, dot product < 0

        all_same = all(vector['x'] == vectors[0]['x'] and vector['y'] == vectors[0]['y'] for vector in vectors)  # Проверка на равенство всех векторов

        if not all_same and prevvector and (vectors[i]['y'] * prevvector['y'] + vectors[i]['x'] * prevvector['x']) < 0:
            # compare magnitude with unit vectors
            vectorlength = math.sqrt(vectors[i]['x']**2 + vectors[i]['y']**2)
            unitv = {'x': vectors[i]['x'] / vectorlength, 'y': vectors[i]['y'] / vectorlength}
            prevlength = math.sqrt(prevvector['x']**2+prevvector['y']**2)
            prevunit = {'x': prevvector['x'] / prevlength, 'y': prevvector['y'] / prevlength}

            # we need to scale down to unit vectors to normalize vector length. Could also just do a tan here
            if abs(unitv['y'] * prevunit['x'] - unitv['x'] * prevunit['y']) < 0.0001:
                continue

        vecd = math.sqrt(vectors[i]['x']**2 + vectors[i]['y']**2)
        if(max_d is None or vecd > max_d):
            max_d = vecd
            translate = vectors[i]
    
    return translate, max_d


def trim_vector(translate, d):
    """
    cuts the vector to a valid moving point
    returns trim vector
    """
    
    vlength2 = translate['x']**2 + translate['y']**2
    if d**2 < vlength2 and not almost_equal(d**2, vlength2):
        scale = math.sqrt((d**2)/vlength2)
        translate['x'] *= scale
        translate['y'] *= scale

    return translate


def find_touching_points(A, B, A_index):
    """
    returns touching points at vertices and on edge segments of polygons
    """

    len_a = len(A['points'])
    len_b = len(B['points'])
    touching = []

    for i in range(len_a-1, -1, -1):
        nexti = len_a-1 if i == 0 else i - 1
        for j in range(len_b):
            nextj = 0 if j == len_b - 1 else j + 1

            if almost_equal(A['points'][i]['x'], B['points'][j]['x'] + B['offsetx']) and almost_equal(
                    A['points'][i]['y'], B['points'][j]['y'] + B['offsety']):
                if(A_index==-1 or A_index == i or (i == len_a-1 and A_index == 0) or (i == A_index-1)):                              
                    A_index = i
                    touching.append({'type': 0, 'A': i, 'B': j})

            elif on_segment(A['points'][i], A['points'][nexti],
                            {'x': B['points'][j]['x']+B['offsetx'], 'y': B['points'][j]['y'] + B['offsety']}):
                touching.append({'type': 1, 'A': nexti, 'B': j})

            elif on_segment(
                    {'x': B['points'][j]['x']+B['offsetx'], 'y': B['points'][j]['y'] + B['offsety']},
                    {'x': B['points'][nextj]['x'] + B['offsetx'], 'y': B['points'][nextj]['y'] + B['offsety']},
                    A['points'][i]):
                touching.append({'type': 2, 'A': i, 'B': nextj})

    return touching, A_index

def is_point_on_line_segment(nfp_point1, nfp_point2, point):
    if not almost_equal((point['x'] - nfp_point1['x']) * (nfp_point2['y'] - nfp_point1['y']), (point['y'] - nfp_point1['y']) * (nfp_point2['x'] - nfp_point1['x'])):
        return False
    
    if min(nfp_point1['x'], nfp_point2['x']) <= point['x'] <= max(nfp_point1['x'], nfp_point2['x']) and min(nfp_point1['y'], nfp_point2['y']) <= point['y'] <= max(nfp_point1['y'], nfp_point2['y']):
        return True
    
    return False

def nfp_polygon(A, B, inside=True):
    """
    given a static polygon A and a movable polygon B, compute a no fit polygon by orbiting B about A
    if the inside flag is set, B is orbited inside of A rather than outside
    if the searchEdges flag is set, all edges of A are explored for NFPs - multiple
    """
    if A is None or len(A['points']) < 3 or B is None or len(B['points']) < 3:
        return None

    # A last point = offsetx, offsety
    A['offsetx'] = 0
    A['offsety'] = 0

    min_a = A['points'][0]['y']
    min_a_index = 0
    max_b = B['points'][0]['y']
    max_b_index = 0

    for i in range(1, len(A['points'])):
        A['points'][i]['marked'] = False
        if A['points'][i]['y'] < min_a:
            min_a = A['points'][i]['y']
            min_a_index = i

    for i in range(1, len(B['points'])):
        B['points'][i]['marked'] = False
        if B['points'][i]['y'] > max_b:
            max_b = B['points'][i]['y']
            max_b_index = i

    if not inside:
        # shift B such that the bottom-most point of B is at the top-most point of A.
        # This guarantees an initial placement with no intersections
        start_point = {
            'x': A['points'][min_a_index]['x'] - B['points'][max_b_index]['x'],
            'y': A['points'][min_a_index]['y'] - B['points'][max_b_index]['y']
        }
    else:
        #  no reliable heuristic for inside
        start_point = search_start_point(A, B, inside)

    NFP_list = []

    while start_point:
        B['offsetx'] = start_point['x']
        B['offsety'] = start_point['y']

        # maintain a list of touching points/edges
        prevvector = None
        NFP = [{
            'x': B['points'][max_b_index]['x'] + B['offsetx'],
            'y': B['points'][max_b_index]['y'] + B['offsety'],
        }]

        reference = {'x': B['points'][max_b_index]['x'] + B['offsetx'],  'y': B['points'][max_b_index]['y'] + B['offsety']}
        startx = reference['x']
        starty = reference['y']
        start_point = None
        
        A_index = -1

        while True:
            touching, A_index = find_touching_points(A, B, A_index)
            vectors = find_feasible_translation_vectors(A, B, touching)
            translate, d = get_single_translation_vector(vectors, prevvector)  

            if translate is None or almost_equal(d, 0):
                if(A_index!=-1):
                    A_index = -1
                    continue
                NFP = None
                break

            translate['start']['marked'] = True
            translate['end']['marked'] = True

            prevvector = translate

            translate = trim_vector(translate, d)

            reference['x'] += translate['x']
            reference['y'] += translate['y']
            
            if almost_equal(reference['x'], startx) and almost_equal(reference['y'], starty):
                # we have made a full loop
                break

            # if A and B start on a touching horizontal line, the end point may not be the start point
            looped = False
            if len(NFP) > 1:
                if is_point_on_line_segment(NFP[-1], reference, NFP[0]):
                    break
                    

            if looped:
                # we've made a full loop
                break

            NFP.append({
                'x': reference['x'],
                'y': reference['y']
            })
            B['offsetx'] += translate['x']
            B['offsety'] += translate['y']

        if NFP and len(NFP) > 0:
            NFP_list.append(NFP)
       
        if(start_point is None):
            start_point = search_start_point(A, B, inside, NFP_list)
        
    return NFP_list

def determine_side(v1, v2):
    """
    A function for determining which side of vector v2 lies vector v1.
    Returns 'right' if v1 is to the right of v2, 'left' if v1 is to the left, and 'on the line',
    if v1 is in a straight line with v2.
    """
    cp = v1[0] * v2[1] - v1[1] * v2[0]
    if cp > 0:
        return 1  # left
    elif cp < 0:
        return -1 # right
    else:
        return 0  # on the line

def choose_translation_vector(vectora, vectorb):
    """
    selects a potential translation vector from two vectors
    """
    vectora = copy.deepcopy(vectora)
    vectora['start'].pop('marked', None)
    vectora['end'].pop('marked', None)

    va_vector = (vectora['end']['x'] - vectora['start']['x'], vectora['end']['y'] - vectora['start']['y'])
    vb_vector = (vectorb['end']['x'] - vectorb['start']['x'], vectorb['end']['y'] - vectorb['start']['y'])

    side = determine_side(vb_vector, va_vector)
    if(vectora['end'] == vectorb['end'] and side == 0):
        return None
    
    if(side == 0):
        return 2

    if(vectora['start'] == vectorb['start']):
        if(side == -1):
            return 1
        if(side == 1):
            return 0
        
    if(vectora['start'] == vectorb['end']):
        if(side == 1):
            return 0
        
    if(vectora['end'] == vectorb['start']):
        if(side == 1):
            return 1
        
    return None

def polygons_intersect_without_edge_touching(poly1, poly2, vector):
    """
    returns true if the polygons intersect but do not touch each other
    """

    d = polygon_slide_distance(poly1, poly2, vector, True)
    vlength2 = vector['x']**2 + vector['y']**2

    if d and d**2 < vlength2 and not almost_equal(d**2, vlength2):
        scale = math.sqrt((d**2)/vlength2)
        vector['x'] *= scale
        vector['y'] *= scale
   
    poly1 = [(point['x'], point['y']) for point in poly1['points']]
    poly2 = [(point['x']+poly2['offsetx']+vector['x'], point['y']+poly2['offsety']+vector['y']) for point in poly2['points']]
    polygon1 = Polygon(poly1)
    polygon2 = Polygon(poly2)
    
    return vector, polygon1.intersects(polygon2) and not polygon1.touches(polygon2)

def search_start_point(A, B, inside=True, NFP=None):
    """
    searches for an arrangement of A and B such that they do not overlap if an NFP is given,
    only search for startpoints that have not already been traversed in the given NFP
    :param A:
    :param B:
    :param inside:
    :param NFP:
    :return:
    """
    # clone arrays
    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    for i in range(0, len(A['points'])-1):
        if not A['points'][i].get('marked'):
            A['points'][i]['marked'] = True
            for j in range(0, len(B['points'])):
                B['offsetx'] = A['points'][i]['x'] - B['points'][j]['x']
                B['offsety'] = A['points'][i]['y'] - B['points'][j]['y']

                bin_side = None
                for k in range(0, len(B['points'])):
                    inpoly = point_in_polygon(
                        {'x': B['points'][k]['x']+B['offsetx'],
                         'y': B['points'][k]['y']+B['offsety']}, A)
                    if inpoly is not None:
                        bin_side = inpoly
                        break

                if bin_side is None:
                    return None

                start_point = {
                    'x': B['offsetx'],
                    'y': B['offsety']
                }
                if ((bin_side and inside) or (not bin_side and not inside)) and (
                            not polygons_intersect_without_edge_touching(A, B, {'x':0, 'y':0})[1] and not inNfp(start_point, NFP)):
                    return start_point

                # slide B along vector
                vx = A['points'][i+1]['x'] - A['points'][i]['x']
                vy = A['points'][i+1]['y'] - A['points'][i]['y']

                d1 = polygon_projection_distance(A, B, {'x': vx, 'y': vy})
                d2 = polygon_projection_distance(B, A, {'x': -vx, 'y': -vy})

                d = None

                if d1 is not None and d2 is not None:
                    d = min(d1, d2)
                elif d1 is None and d2 is not None:
                    d = d2
                elif d1 is not None and d2 is None:
                    d = d1

                # only slide until no longer negative
                if not (d is not None and not almost_equal(d, 0) and d > 0):
                    continue

                vd2 = vx * vx + vy * vy
                if d * d < vd2 and not almost_equal(d*d, vd2):
                    vd = math.sqrt(vx * vx + vy * vy)
                    vx *= d /vd
                    vy *= d /vd

                B['offsetx'] += vx
                B['offsety'] += vy

                for k in range(0, len(B['points'])):
                    inpoly = point_in_polygon(
                        {'x': B['points'][k]['x']+B['offsetx'],
                         'y': B['points'][k]['y']+B['offsety']}, A)
                    if inpoly is not None:
                        bin_side = inpoly
                        break

                start_point = {'x': B['offsetx'], 'y': B['offsety']}
                if ((bin_side and inside) or (not bin_side and not inside)) and (
                            not polygons_intersect_without_edge_touching(A, B, {'x':0, 'y':0})[1] and not inNfp(start_point, NFP)):
                    return start_point

    return None


def inNfp(p, nfp):
    """
    returns true if point already exists in the given nfp
    :param p:
    :param nfp:
    :return:
    """
    if not nfp or len(nfp) == 0:
        return False

    for i in range(0, len(nfp)):
        for j in range(0, len(nfp[i])):
            if almost_equal(p['x'], nfp[i][j]['x']) and almost_equal(p['y'], nfp[i][j]['y']):
                return True

    return False


def point_in_polygon(point, polygon):
    if isinstance(polygon, list):
        polygon = {'points': polygon}
    if len(polygon.get('points')) < 3:
        return None

    inside = False
    offsetx = polygon.get('offsetx') or 0
    offsety = polygon.get('offsety') or 0

    j = len(polygon['points']) - 1
    for i in range(0, len(polygon['points'])):
        xi = polygon['points'][i]['x'] + offsetx
        yi = polygon['points'][i]['y'] + offsety
        xj = polygon['points'][j]['x'] + offsetx
        yj = polygon['points'][j]['y'] + offsety

        if almost_equal(xi, point['x']) and almost_equal(yi, point['y']):
            return None

        if on_segment({'x': xi, 'y': yi}, {'x':xj, 'y':yj}, point):
            return None  # exactly on the segment

        if almost_equal(xi, xj) and almost_equal(yi, yj):
            # ignore very small lines
            continue

        intersect = ((yi > point['y']) != (yj > point['y'])) and (point['x'] < (xj - xi) * (point['y'] - yi) / (yj - yi) + xi)
        if intersect:
            inside = not inside

    return inside


def polygon_projection_distance(A, B, direction):
    """
    project each point of B onto A in the given direction, and return the distance
    :param A:
    :param B:
    :param direction:
    :return:
    """
    b_offsetx = B.get('offsetx') or 0
    b_offsety = B.get('offsety') or 0
    a_offsetx = A.get('offsetx') or 0
    a_offsety = A.get('offsety') or 0

    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    edge_a = A['points']
    edge_b = B['points']
    distance = None
    p = dict()
    s1 = dict()
    s2 = dict()
    for i in range(0, len(edge_b)):
        # the shortest/most negative projection of B onto A
        min_projection = minp = None
        for j in range(0, len(edge_a) - 1):
            p['x'] = edge_b[i]['x'] + b_offsetx
            p['y'] = edge_b[i]['y'] + b_offsety
            s1['x'] = edge_a[j]['x'] + a_offsetx
            s1['y'] = edge_a[j]['y'] + a_offsety
            s2['x'] = edge_a[j+1]['x'] + a_offsetx
            s2['y'] = edge_a[j+1]['y'] + a_offsety

            if abs((s2['y'] - s1['y']) * direction['x'] - (s2['x'] - s2['x']) * direction['y']) < TOL_ERROR:
                continue

            # project point, ignore edge boundaries
            d = point_distance(p, s1, s2, direction)
            if d and (min_projection is None or d < min_projection):
                min_projection = d

        if min_projection and (distance is None or min_projection > distance):
            distance = min_projection

    return distance


def point_distance(p, s1, s2, normal, infinite=None):
    normal = normalize_vector(normal)

    dir_point = {
        'x': normal['y'],
        'y': -normal['x'],
    }

    pdot = p['x'] * dir_point['x'] + p['y'] * dir_point['y']
    s1dot = s1['x'] * dir_point['x'] + s1['y'] * dir_point['y']
    s2dot = s2['x'] * dir_point['x'] + s2['y'] * dir_point['y']

    pdotnorm = p['x']*normal['x'] + p['y'] * normal['y']
    s1dotnorm = s1['x']*normal['x'] + s1['y'] * normal['y']
    s2dotnorm = s2['x'] * normal['x'] + s2['y'] * normal['y']

    if infinite is None:
        # dot doesn't collide with segment, or lies directly on the vertex
        if ((pdot < s1dot or almost_equal(pdot, s1dot)) and (pdot < s2dot or almost_equal(pdot, s2dot))) or (
                    (pdot > s1dot or almost_equal(pdot, s1dot)) and ((pdot > s2dot) or almost_equal(pdot, s2dot))):
            return None
        if (almost_equal(pdot, s1dot) and almost_equal(pdot, s2dot)) and (
                        pdotnorm > s1dotnorm and pdotnorm > s2dotnorm):
            return min(pdotnorm - s1dotnorm, pdotnorm - s2dotnorm)
        if almost_equal(pdot, s1dot) and almost_equal(pdot, s2dot) and pdotnorm < s1dotnorm and pdotnorm < s2dotnorm:
            return -min(s1dotnorm-pdotnorm, s2dotnorm-pdotnorm)

    return -(pdotnorm - s1dotnorm + (s1dotnorm - s2dotnorm) * (s1dot - pdot)/(s1dot - s2dot))


def polygon_slide_distance(A, B, direction, ignorenegative):

    b_offsetx = B.get('offsetx') or 0
    b_offsety = B.get('offsety') or 0
    a_offsetx = A.get('offsetx') or 0
    a_offsety = A.get('offsety') or 0

    A = copy.deepcopy(A)
    B = copy.deepcopy(B)

    if not A['points'][-1] == A['points'][0]:
        A['points'].append(A['points'][0])
    if not B['points'][0] == B['points'][-1]:
        B['points'].append(B['points'][0])

    edge_a = A['points']
    edge_b = B['points']
    distance = None

    dir_point = normalize_vector(direction)

    for i in range(0, len(edge_b) - 1):

        for j in range(0, len(edge_a) - 1):
            A1 = {'x': edge_a[j]['x'] + a_offsetx, 'y': edge_a[j]['y'] + a_offsety}
            A2 = {'x': edge_a[j+1]['x'] + a_offsetx, 'y': edge_a[j+1]['y'] + a_offsety}
            B1 = {'x': edge_b[i]['x'] + b_offsetx, 'y': edge_b[i]['y'] + b_offsety}
            B2 = {'x': edge_b[i + 1]['x'] + b_offsetx, 'y': edge_b[i + 1]['y'] + b_offsety}

            if (almost_equal(A1['x'], A2['x']) and almost_equal(A1['y'], A2['y'])) or almost_equal(
                    B1['x'], B2['x']) and almost_equal(B1['y'], B2['y']):
                continue

            d = segment_distance(A1, A2, B1, B2, dir_point)
            if d and (distance is None or d < distance):
                if not ignorenegative or d > 0 or almost_equal(d, 0):
                    distance = d
    return distance


def segment_distance(A, B, E, F, direction):
    normal = {
        'x': direction['y'],
        'y': -direction['x']
    }
    reverse = {
        'x': -direction['x'],
        'y': -direction['y']
    }

    dot_a = A['x'] * normal['x'] + A['y'] * normal['y']
    dot_b = B['x'] * normal['x'] + B['y'] * normal['y']
    dot_e = E['x'] * normal['x'] + E['y'] * normal['y']
    dot_f = F['x'] * normal['x'] + F['y'] * normal['y']

    cross_a = A['x'] * direction['x'] + A['y'] * direction['y']
    cross_b = B['x'] * direction['x'] + B['y'] * direction['y']
    cross_e = E['x'] * direction['x'] + E['y'] * direction['y']
    cross_f = F['x'] * direction['x'] + F['y'] * direction['y']

    ab_min = min(dot_a, dot_b)
    ab_max = max(dot_a, dot_b)

    ef_min = min(dot_e, dot_f)
    ef_max = max(dot_e, dot_f)

    # segments that will touch at one point
    if almost_equal(ab_max, ef_min, TOL_ERROR) or almost_equal(ab_min, ef_max, TOL_ERROR):
        return None

    # segments miss each other completely
    if ab_max < ef_min or ab_min > ef_max:
        return None

    if (ab_max > ef_max and ab_min < ef_min) or (ef_max > ab_max and ef_min < ab_min):
        overlap = 1
    else:
        min_max = min(ab_max, ef_max)
        max_min = max(ab_min, ef_min)
        max_max = max(ab_max, ef_max)
        min_min = min(ab_min, ef_min)

        overlap = (min_max - max_min) / (max_max - min_min)

    cross_abe = (E['y'] - A['y']) * (B['x'] - A['x']) - (E['x'] - A['x']) * (B['y'] - A['y'])
    cross_abf = (F['y'] - A['y']) * (B['x'] - A['x']) - (F['x'] - A['x']) * (B['y'] - A['y'])

    # lines are colinear
    if almost_equal(cross_abe, 0) and almost_equal(cross_abf, 0):
        ab_norm = {'x': B['y'] - A['y'], 'y': A['x'] - B['x']}
        ef_norm = {'x': F['y'] - E['y'], 'y': E['x'] - F['x']}

        ab_norm_length = math.sqrt(ab_norm['x']**2 + ab_norm['y']**2)
        ab_norm['x'] /= ab_norm_length
        ab_norm['y'] /= ab_norm_length

        ef_norm_length = math.sqrt(ef_norm['x']**2 + ef_norm['y']**2)
        ef_norm['x'] /= ef_norm_length
        ef_norm['y'] /= ef_norm_length

        # segment normals must point in opposite directions
        if abs(ab_norm['y'] * ef_norm['x'] - ab_norm['x'] * ef_norm['y']) < TOL_ERROR and (
                            ab_norm['y'] * ef_norm['y'] + ab_norm['x'] * ef_norm['x'] < 0):
            # normal of AB segment must point in same direction as given direction vector
            norm_dot = ab_norm['y'] * direction['y'] + ab_norm['x'] * direction['x']
            # the segments merely slide along eachother
            if almost_equal(norm_dot, 0, TOL_ERROR):
                return None

            if norm_dot < 0:
                return 0

        return None

    distances = list()

    # coincident points
    if almost_equal(dot_a, dot_e):
        distances.append(cross_a - cross_e)
    elif almost_equal(dot_a, dot_f):
        distances.append(cross_a - cross_f)
    elif ef_min < dot_a and dot_a < ef_max:
        d = point_distance(A, E, F, reverse)
        # A currently touches EF, but AB is moving away from EF
        if d and almost_equal(d, 0):
            db = point_distance(B, E, F, reverse, True)
            if db < 0 or almost_equal(db * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if almost_equal(dot_b, dot_e):
        distances.append(cross_b - cross_e)
    elif almost_equal(dot_b, dot_f):
        distances.append(cross_b - cross_f)
    elif dot_b > ef_min and dot_b < ef_max:
        d = point_distance(B, E, F, reverse)

        if d and almost_equal(d, 0):
            da = point_distance(A, E, F, reverse, True)
            if da < 0 or almost_equal(da * overlap, 0):
                d = None

        if d:
            distances.append(d)

    if dot_e > ab_min and dot_e < ab_max:
        d = point_distance(E, A, B, direction)
        if d and almost_equal(d, 0):
            df = point_distance(F, A, B, direction, True)
            if df < 0 or almost_equal(df * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if len(distances) == 0:
        return None

    if dot_f > ab_min and dot_f < ab_max:
        d = point_distance(F, A, B, direction)
        if d and almost_equal(d, 0):
            de = point_distance(E, A, B, direction, True)
            if de < 0 or almost_equal(de * overlap, 0):
                d = None
        if d:
            distances.append(d)

    if len(distances) == 0:
        return None

    return min(distances)


def polygon_area(polygon):
    area = 0
    j = len(polygon) - 1
    for i in range(0, len(polygon)):
        area += (polygon[j]['x'] + polygon[i]['x']) * (polygon[j]['y'] - polygon[i]['y'])
        j = i

    return 0.5 * area



if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as matplotlibPolygon

    # Тест 1
    A = {'points':[{'x': 6, 'y':2}, {'x': 8, 'y':4}, {'x': 10, 'y':2}, {'x': 8, 'y':0}]}
    B = {'points':[{'x': 0, 'y':0}, {'x':2, 'y':2}, {'x':4, 'y':0} ]}

    # Test 2
    A = {'points':[{'x': 5, 'y':7}, {'x': 8, 'y':7}, {'x': 8, 'y':5}, {'x': 6, 'y':5},
                   {'x': 6, 'y':3}, {'x': 8, 'y':3}, {'x': 8, 'y':1}, {'x': 10, 'y':1},
                   {'x': 10, 'y':3}, {'x': 12, 'y':3}, {'x': 12, 'y':5}, {'x': 10, 'y':5},
                   {'x': 10, 'y':7}, {'x': 13, 'y':7}, {'x': 13, 'y':0}, {'x': 5, 'y':0}]}
    B = {'points':[{'x': 0, 'y':2}, {'x':2, 'y':2}, {'x':2, 'y':0}, {'x':0, 'y':0}]}

    # Test 3
    A = {'points':[{'x': 0, 'y':0}, {'x': 0, 'y':6}, {'x': 5, 'y':5.5}, {'x': 7, 'y':3.5},
                   {'x': 3, 'y':5}, {'x': 1, 'y':2.5}, {'x': 4, 'y':1}, {'x': 6, 'y':2},
                   {'x': 6, 'y':0}, ]}
    B = {'points':[{'x': -6, 'y':-3}, {'x':-5, 'y':-1}, {'x':-2, 'y':-1}, {'x':-1, 'y':-3.5},
                   {'x': -2, 'y':-2}, {'x':-4, 'y':-1.5}, ]}
    
    # Test 4
    A = {'points': [{'x': 6, 'y': 2}, {'x': 6, 'y': 4}, {'x': 8, 'y': 4}, {'x': 8, 'y': 2}]}
    B = {'points': [{'x': 0, 'y': 0}, {'x': 0, 'y': 2}, {'x': 2, 'y': 2}, {'x': 2, 'y': 0}]}

    # Test 5
    A = {'points': [{'x': 0.0, 'y': 2.0}, {'x': 7.0, 'y': 0.0}, {'x': 8.0, 'y': 2.0}, {'x': 8.0, 'y': 3.0},{'x': 9.0, 'y': 4.0}, {'x': 8.0, 'y': 5.0}, {'x': 7.0, 'y': 7.0}, {'x': 0.0, 'y': 6.0}]}
    B = {'points': [{'x': 0.0, 'y': -7.0}, {'x': 1.0, 'y': -8.0}, {'x': 2.0, 'y': -10.0}, {'x': 9.0, 'y': -9.0}, {'x': 9.0, 'y': -5.0}, {'x': 2.0, 'y': -3.0}, {'x': 1.0, 'y': -5.0}, {'x': 1.0, 'y': -6.0}]}

    if polygon_area(A['points']) < 0:
        A['points'].reverse()
    if polygon_area(B['points']) < 0:
        B['points'].reverse()

    fig, ax = plt.subplots()
    polygonA = matplotlibPolygon([(path['x'], path['y']) for path in A['points']], closed=True, fill='#BED574', facecolor='#BED574', edgecolor='#72B73E')
    polygonB = matplotlibPolygon([(path['x'], path['y']) for path in B['points']], closed=True, fill='#E8B0A3', facecolor='#E8B0A3', edgecolor='#EA562B')
    ax.add_patch(polygonA)
    ax.add_patch(polygonB)

    # for outer NFPs, the first is guaranteed to be the largest.
    # Any subsequent NFPs that lie inside the first are hole
    nfp = nfp_polygon(A, B, False)

    for i in range(0, len(nfp)):
        if polygon_area(nfp[i]) > 0:
            nfp[i].reverse()

        if i > 0:
            if point_in_polygon(nfp[i][0], nfp[0]):
                if polygon_area(nfp[i]) < 0:
                    nfp[i].reverse()

    x_list_A = [path['x'] for path in A['points']]
    y_list_A = [path['y'] for path in A['points']]
    x_list_B = [path['x'] for path in B['points']]
    y_list_B = [path['y'] for path in B['points']]
    max_x = max(max(x_list_A),max(x_list_B))
    min_x = min(min(x_list_A),min(x_list_B))
    max_y = max(max(y_list_A),max(y_list_B))
    min_y = min(min(y_list_A),min(y_list_B))

    for nf in nfp:
        path_list = []

        for path in nf:
            path_list.append((path['x'], path['y']))
            if(min_x > path['x']):
                min_x = path['x']
            elif(max_x < path['x']):
                max_x = path['x']
            if(min_y > path['y']):
                min_y = path['y']
            elif(max_y < path['y']):
                max_y = path['y']
       
        polygon = matplotlibPolygon(path_list, closed=True, fill=None, edgecolor='#19692F')
        polygon.set_linewidth(2)

        ax.add_patch(polygon)

    plt.xlim(min_x-1, max_x+1)
    plt.ylim(min_y-1, max_y+1)

    plt.show()