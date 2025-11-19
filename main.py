import warnings
import osmnx as ox # v 1.9.4
import networkx as nx # v 3.3
import numpy as np # v 1.26.4
import matplotlib.pyplot as plt # v 3.9.2
import geopandas as gpd # v 0.14.4
from shapely.geometry import MultiLineString # v 2.0.6
import math
import scipy.io as sio
#import algorithms

# ******************************************************
warnings.filterwarnings("ignore", category=FutureWarning,
                        message="The expected order of coordinates in `bbox` will change in the v2.0.0 release")

# ******************************************************

def graph_extraction(bbox):
    # create a graph with given bbox
    G_initial = ox.graph_from_bbox(bbox=bbox, # select area of interest
                            network_type='drive', # filter only drive streets
                            simplify=False, # simplify topology
                            retain_all=False,  # retain only the largest weakly connected component
                            truncate_by_edge=True)  # retain nodes outside bbox if at least one of node’s neighbors is within the bounding box
    
    # extract nodes (intersections) and edges (roads) from the graph
    nodes_G_initial, edges_G_initial = ox.graph_to_gdfs(G_initial) # convert a MultiDiGraph to node and edge GeoDataFrames

    print(f"G_initial: {G_initial}")  # MultiDiGraph
    print(f"Number of nodes in G_initial: {len(nodes_G_initial)}")  
    print(f"NUmber of edges in G_initial: {len(edges_G_initial)}")  

    return G_initial 

def graph_processing(G):
    
    # SIMPLIFY GRAPH REMOVING INTERSTITIAL NODES
    G_simpl = ox.simplify_graph(G)

    # PROJECT THE GRAPH ON UTM COORDINATES
    G_simpl_proj = ox.project_graph(G_simpl)

    edge_centrality = nx.closeness_centrality(nx.line_graph(G_simpl_proj))

    nx.set_edge_attributes(G_simpl_proj, edge_centrality, 'edge_centrality')

    nodes_G_simpl_proj, edges_G_simpl_proj = ox.graph_to_gdfs(G_simpl_proj, nodes= True, edges = True) # convert a MultiDiGraph to node and edge GeoDataFrames (geopandas)

    print(f"Number of nodes in G_simpl_proj: {len(nodes_G_simpl_proj)}")  
    print(f"Number of edges in G_simpl_proj: {len(edges_G_simpl_proj)}") 

    # FILTER OUT CUL-DE-SACS
    # get DegreeView object of the graph: iterator for (node, degree)
    degrees = G_simpl_proj.degree()
    nodes_to_remove = [] 
    for index, row in edges_G_simpl_proj.iterrows():
        # consider only streets (edges) of low importance
        if edges_G_simpl_proj.loc[index,'highway'] != 'primary' or edges_G_simpl_proj.loc[index,'highway'] != 'secondary':
            # extract the two extremes of the edge
            u = index[0]
            v = index[1]
            # check if one of the two extremes is not connected to any other edge 
            if degrees[u] == 1:
                nodes_to_remove.append(u)
            if degrees[v] == 1: 
                nodes_to_remove.append(v) 
            if degrees[u] == 2:
                in_edges_u = []
                for e_in_u in G_simpl_proj.in_edges(u):
                    in_edges_u.append(e_in_u)
                out_edges_u=[]
                for e_out_u in G_simpl_proj.out_edges(u):
                    out_edges_u.append(e_out_u)
                if in_edges_u and out_edges_u:
                    if set(in_edges_u[0]) == set(out_edges_u[0]):
                        nodes_to_remove.append(u)
                        # edges_to_remove[index] = 1
            if degrees[v] == 2:
                in_edges_v = []
                for e_in_v in G_simpl_proj.in_edges(v):
                    in_edges_v.append(e_in_v)
                out_edges_v=[]
                for e_out_v in G_simpl_proj.out_edges(v):
                    out_edges_v.append(e_out_v)
                if in_edges_v and out_edges_v:
                    if set(in_edges_v[0]) == set(out_edges_v[0]):
                        nodes_to_remove.append(v)
                        # edges_to_remove[index] = 1
                
    G_filt = G_simpl_proj.copy()
    G_filt.remove_nodes_from(nodes_to_remove)
    return G_filt


def get_adj(G):
    # get graph adjacency matrix
    Adj = nx.to_numpy_array(G)
    return Adj

def get_borders(edges):
    # get if each edge is border or not
    borders = []
    for index, row in edges.iterrows():
        if edges.loc[index, 'isBorder'] == True:
            borders.append(edges.loc[index, 'highway'])
        else:
            borders.append(0)
    return borders

def get_categories(edges):
    # get the category associated to each edge
    categories = []
    for index, row in edges.iterrows():
        categories.append(edges.loc[index, 'highway'])
    return categories


def graph_visualization(G, title = "", want_to_save = False):
    # visualize ...
    fig, ax = ox.plot_graph(G, bgcolor='k', node_size=30, node_color='#999999',
                            node_edgecolor='none', node_zorder=2,
                            edge_color='#555555', edge_linewidth=2, edge_alpha=1,
                            show=False, close=False)
    fig.suptitle(title, color='k')
    if want_to_save:
        fig.savefig('graph\Original Graph.png', dpi=500, bbox_inches='tight')
    plt.show()

def remove_duplicates(l):
    #     INPUT:
    #         - l : list
    #     OUTPUT:
    #         - l_without_duplicates : l without any duplicate element

    l_without_duplicates = []
    for e in l:
        if e not in l_without_duplicates:
            l_without_duplicates.append(e)
    return l_without_duplicates

def nodes_intersections_list(geoDf):
    #     INPUT:
    #         - geoDf : GeoDataFrame in which the, at each index, the labels 'u'
    #                   and 'v' contain lists of source/destination nodes
    #     OUTPUT:
    #         - intersections_list : list such that the i-th entry contains the
    #                                list of all nodes ('u' or 'v') of i-th index
    #                                with at least a correspondent in any other
    #                                index of geoDf

    intersections_list = []
    for index, row in geoDf.iterrows():

        # extract all nodes in current index
        nodes_actual = geoDf.loc[index, 'u'] + geoDf.loc[index, 'v']
        nodes_actual = remove_duplicates(nodes_actual)

        # create a list with all nodes of all the other indeces
        nodes_others = []
        for i, r in geoDf.iterrows():
            if i != index:
                nodes_others = nodes_others + geoDf.loc[i, 'u'] + geoDf.loc[i, 'v']

        nodes_others = remove_duplicates(nodes_others)

        # create a list of all intersection points between current street
        # and the other street
        intersections = [node \
                         for node in nodes_actual \
                         if node in nodes_others]

        intersections_list.append(intersections)

    return intersections_list

def category_score_map(category):
    #     INPUT:
    #         - category : cateogry of the road
    #     OUTPUT:
    #         - score : score associated to the given category

    all_possible_categories = ['motorway', 'motorway_link', \
                               'trunk', 'trunk_link', \
                               'primary', 'primary_link', \
                               'secondary', 'secondary_link', \
                               'tertiary', 'tertiary_link', \
                               'residential', 'unclassified']
    # there are many others more specific

    number_of_classes = 7
    all_possibile_scores = list(range(number_of_classes)[::-1])  # [6, 5, 4, 3, 2, 1, 0]

    score = -1
    for i in range(len(all_possible_categories)):
        if any(category in s for s in all_possible_categories[i:i + 1]):
            score = all_possibile_scores[math.trunc(i / 2)]
    if score == -1:
        score = all_possibile_scores[-1]

    return score

def create_LTI_system(Adj, categories, borders):
    '''
    Consider the following LTI system
    x(t+1) = Ax(t) + Bu(t)
    y(t) = Cx(t)
    '''
    n = np.shape(Adj)[0] # number of roads in the network --> system dimension
    p = 8 # number of sensors to choose 
    C = np.eye(n) # all possible measurements 
    A= Adj + np.eye(n)
    B = np.zeros((n,1))
    for i in range(n):
        v = np.zeros(n)
        for j in range(n):
            if Adj[i, j] == 1 or i == j:
                v[j] = categories[j]
        
        if borders[i] == 0:
            A[np.arange(n) != i, i] = v[np.arange(n) != i] / np.sum(v)
            A[i, i] = 1 - np.sum(A[np.arange(n) != i, i])
        else:
            A[np.arange(n) != i, i] = v[np.arange(n) != i] / (np.sum(v) + v[i])
            A[i, i] = 1 - np.sum(A[np.arange(n) != i, i]) - v[i] / (np.sum(v) + v[i])
            B[i] = v[i] / (np.sum(v) + v[i])
    return A,B,C



# ******************************************************


if __name__ == "__main__":
    # bounding box latitude-longitude coordinates to select the area of interest
    bbox = [45.41080,45.40653,11.90678,11.89139] # stanga
    # bbox = [45.43306, 45.42866, 11.89090, 11.88072]  # arcella
    # bbox = [45.43306, 45.42815, 11.89507, 11.88072]  # arcella (more long distance)
    # bbox = [45.40000,45.38418,11.95493,11.92853] # padova industrial zone
    G_init = graph_extraction(bbox)
    graph_visualization(G_init, title='Original Graph (simplified version)')
    G_processed = graph_processing(G_init)
    # extract nodes and edges 
    nodes_processed, edges_processed = ox.graph_to_gdfs(G_processed, nodes=True, edges=True)
    # identify which edge is on the border of the area of interest
    borders = []
    degree_dict = G_processed.degree()
    for index, row in edges_processed.iterrows():
        # extract the two extremes of the actual edge
        u = index[0]
        v = index[1]
        # check if one of the two extremes is not connected to any other edge
        f = 0
        if degree_dict[u] == 1:
            f = 1
        elif degree_dict[v] == 1:
            f = 1
        elif degree_dict[u] == 2:
            in_edges_u = []
            for e_in_u in G_processed.in_edges(u):
                in_edges_u.append(e_in_u)
            out_edges_u = []
            for e_out_u in G_processed.out_edges(u):
                out_edges_u.append(e_out_u)
            if in_edges_u and out_edges_u:
                if set(in_edges_u[0]) == set(out_edges_u[0]):
                    f = 1
        elif degree_dict[v] == 2:
            in_edges_v = []
            for e_in_v in G_processed.in_edges(v):
                in_edges_v.append(e_in_v)
            out_edges_v = []
            for e_out_v in G_processed.out_edges(v):
                out_edges_v.append(e_out_v)
            if in_edges_v and out_edges_v:
                if set(in_edges_v[0]) == set(out_edges_v[0]):
                    f = 1
        if f == 0:
            borders.append(False)
        else:
            borders.append(True)
    edges_processed['isBorder'] = borders

    # remove duplicate streets
    edges_preproc = edges_processed.copy()
    for index, row in edges_processed.iterrows():
        # consider the edge characterized by nodes (u,v) and given osmid
        u = index[0]
        v = index[1]
        osmid = edges_processed.loc[index, 'osmid']
        # scan all the next ones
        for i, r in edges_processed.iterrows():
            if i > index:
                # if one row has the same (u,v) or (v,u)...
                if (i[0] == u and i[1] == v) or (i[0] == v and i[1] == u):
                    # ... and same osmID...
                    if edges_processed.loc[i, 'osmid'] == osmid:
                        # ... remove it (it is a duplicate)
                        edges_preproc = edges_preproc.drop([i], axis=0)
    # there are somea edges with multiple nmes -> consider only the first one
    edges_preproc = edges_preproc.sort_index()
    for index, row in edges_preproc.iterrows():
        if type(row['name']) == list:
            edges_preproc.loc[index, 'name'] = row['name'][0]

    # group the edges by street name.  
    name_groups = edges_preproc.groupby('name')
    edges_grouped = gpd.GeoDataFrame(data=None, index=range(name_groups.ngroups),
                                 columns=['osmid', 'oneway', 'lanes', 'name', 'highway', 'reversed', 'length',
                                          'geometry', 'edge_centrality', 'maxspeed', 'junction', 'isBorder', 'u', 'v'])

    multiline = []
    i = 0
    for name, group in name_groups:
        # name: name of the street
        osmid_list = []
        oneway_list = []
        lanes_list = []
        name_list = []
        highway_list = []
        reversed_list = []
        length_list = []
        geometry_list = []
        edge_centrality_list = []
        maxspeed_list = []
        junction_list = []
        isBorder_list = []
        u_list = []
        v_list = []
        for index, row in group.iterrows():  
            u_list.append(index[0])  # We get value from column 'u'
            v_list.append(index[1])  # We get value from column 'v'
            osmid_list.append(row['osmid'])
            oneway_list.append(row['oneway'])
            lanes_list.append(row['lanes'])
            highway_list.append(row['highway'])
            reversed_list.append(row['reversed'])
            length_list.append(row['length'])
            geometry_list.append(row['geometry'])
            edge_centrality_list.append(row['edge_centrality'])
            maxspeed_list.append(row['maxspeed'])
            junction_list.append(row['junction'])
            isBorder_list.append(row['isBorder'])
        edges_grouped.at[i, 'u'] = u_list
        edges_grouped.at[i, 'v'] = v_list
        edges_grouped.at[i, 'osmid'] = osmid_list
        edges_grouped.at[i, 'oneway'] = oneway_list
        edges_grouped.at[i, 'lanes'] = lanes_list
        edges_grouped.at[i, 'name'] = name
        edges_grouped.at[i, 'highway'] = highway_list
        edges_grouped.at[i, 'reversed'] = reversed_list
        edges_grouped.at[i, 'length'] = length_list
        # edges_grouped.at[i, 'geometry']=geometry_list
        multiline.append(MultiLineString(geometry_list))
        edges_grouped.at[i, 'edge_centrality'] = edge_centrality_list
        edges_grouped.at[i, 'maxspeed'] = maxspeed_list
        edges_grouped.at[i, 'junction'] = junction_list
        edges_grouped.at[i, 'isBorder'] = isBorder_list
        i += 1
    edges_grouped['geometry'] = multiline

    intersections_list = nodes_intersections_list(edges_grouped)

    i = 0
    for index, row in edges_grouped.iterrows():
        cols = list(edges_grouped)
        for label in cols:
            if label == 'edge_centrality':
                edges_grouped.at[index, label] = \
                    np.mean(np.asarray(edges_grouped.at[index, label]))
            elif label == 'highway':
                l_strings = edges_grouped.at[index, label]
                l_scores = []
                for j in range(len(l_strings)):
                    if type(l_strings[j]) is list:
                        for k in range(len(l_strings[j])):
                            l_scores.append(category_score_map(l_strings[j][k]))
                    else:
                        l_scores.append(category_score_map(l_strings[j]))
                edges_grouped.at[index, label] = max(l_scores)
            elif label == 'geometry':
                edges_grouped.at[index, label] = MultiLineString(edges_grouped.at[index, label])
            elif label == 'length':
                edges_grouped.at[index, label] = \
                    np.sum(np.asarray(edges_grouped.at[index, label]))
            elif label == 'isBorder':
                if True in edges_grouped.at[index, label]:
                    edges_grouped.at[index, label] = True
                else:
                    edges_grouped.at[index, label] = False
            elif label == 'name':
                if type(edges_grouped.at[index, label]) == list:
                    edges_grouped.at[index, label] = edges_grouped.at[index, label][0]
            elif label == 'osmid':
                # just eliminate duplicates
                edges_grouped.at[index, label] = \
                    remove_duplicates(edges_grouped.at[index, label])
            elif label == 'u' or label == 'v':
                # just eliminate duplicates
                edges_grouped.at[index, label] = \
                    remove_duplicates(edges_grouped.at[index, label])
                # and remove all internal nodes (that does not intersect with any
                # other street)
                edges_grouped.at[index, label] = \
                    set(edges_grouped.at[index, label]) & \
                    set(intersections_list[i])
                edges_grouped.at[index, label] = list(edges_grouped.at[index, label])
        i = i + 1

    G = nx.Graph()
    for index, row in edges_grouped.iterrows():
        # add a node with the current street
        G.add_node(index, name=edges_grouped.loc[index, 'name'])
        nodes_actual = edges_grouped.loc[index, 'u'] + edges_grouped.loc[index, 'v']
        nodes_actual = remove_duplicates(nodes_actual)
        # explore all the other streets
        for i, r in edges_grouped.iterrows():
            if i != index:
                # create a list of all nodes of another street
                nodes_other = edges_grouped.loc[i, 'u'] + edges_grouped.loc[i, 'v']
                nodes_other = remove_duplicates(nodes_other)

                # create a list of all intersection points between current street
                # and the other street
                intersections = [node \
                                for node in nodes_actual \
                                if node in nodes_other]
                # if there is at least one intersection...
                if len(intersections) >= 1:
                    # ... create an edge between the two streets (according to the
                    # dual graph convention)
                    G.add_edge(index, i, crossRoad='' + \
                                                edges_grouped.loc[index, 'name'] + \
                                                ' - ' + edges_grouped.loc[i, 'name'] + '')

    degree_G = G.degree()
    names_to_remove = []
    for n, d in degree_G:
        if d == 1:
            names_to_remove.append(G.nodes[n]['name'])
    edges_final = edges_grouped.copy()

    Adj = get_adj(G)
    borders = get_borders(edges_final)
    categories = get_categories(edges_final)
    A,B,C = create_LTI_system(Adj, categories, borders)

    #np.savetxt('.\A.mat', A)
    #np.savetxt('.\B.mat', B)
    #np.savetxt('.\C.mat', C)
    sio.savemat('.\A.mat', {'A': A})
    sio.savemat('.\B.mat', {'B': B})
    sio.savemat('.\C.mat', {'C': C})