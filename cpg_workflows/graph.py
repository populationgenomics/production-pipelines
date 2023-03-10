"""
Code to display the prod pipelines graph

Heavily inspired by:
https://towardsdatascience.com/visualize-hierarchical-data-using-plotly-and-datapane-7e5abe2686e1
"""

import os
from copy import deepcopy
from typing import Callable
from itertools import groupby

import networkx as nx
from plotly.subplots import make_subplots
import plotly.graph_objs as go


class GraphPlot:
    def __init__(self, G: nx.DiGraph, **kwargs):
        # Data
        self.G = deepcopy(G)

        script_location = os.path.dirname(os.path.abspath(__file__))
        self.stages_folder = os.path.abspath(f'{script_location}/stages/')
        self.links = {n: self._get_node_file_link(n) for n in G.nodes}

        # Text and sizes
        self.title = 'Workflow Graph'
        self.title_fontsize = 24
        self.node_text_position = 'top center'
        self.node_text_fontsize = 14
        self.node_size = 50
        self.node_border_weight = 5
        self.edge_weight = 1

        # Layout
        self.partite_key = 'layer'
        self.align = 'horizontal'
        self.layout_scale = 10
        self.show_legend = False

        # Colors
        self.colorscale = 'Blugrn'
        self.node_color_key = self.partite_key
        self.grey_color = 'rgba(153, 153, 153, 0.5)'
        self.dark_emphasis_color = '#0c1f27'
        self.arrow_opacity = 0.8

        # Convert any kwargs into attributes
        self.__dict__.update(kwargs)

        # Graph transforms
        # Add name as attribute
        nx.set_node_attributes(self.G, {n: n for n in self.G.nodes}, name='name')

        # Reverse all edges
        self.G = G.reverse()

        # Recalculate the depths using the topological order
        self._recalculate_depth(new_key=self.partite_key)

        # Calculate the depth_order and position
        self._calculate_depth_order(layer_key=self.partite_key, new_key='layer_order')

    def __add__(self, other) -> go.Figure:
        assert type(self) == type(other)

        fig = make_subplots(rows=1, cols=2, column_titles=[self.title, other.title])

        for trace in self.create_traces():
            fig.add_trace(trace, row=1, col=1)

        for trace in other.create_traces():
            fig.add_trace(trace, row=1, col=2)

        # Set the overall layout using this objects layout
        # Make sure the updates are to both subplot axes
        layout = self._get_layout()
        fig.update_layout(layout)
        fig.update_layout(title='')
        fig.update_yaxes(layout.yaxis)
        fig.update_xaxes(layout.xaxis)

        fig.layout.update(
            annotations=self.get_annotations(xref='x1', yref='y1')
            + other.get_annotations(xref='x2', yref='y2')
        )

        return fig

    def _get_layout(self):
        return go.Layout(
            title=self.title,
            titlefont_size=self.title_fontsize,
            showlegend=self.show_legend,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(
                showgrid=False,
                zeroline=False,
                showticklabels=False,
                autorange='reversed',
            ),
            coloraxis=dict(reversescale=False),
        )

    def display_graph(self):
        fig = self.create_figure()
        fig.show()

    def create_traces(self) -> go.Figure:
        # Add weight and depth attributes to the nodes
        for node in self.G.nodes:
            self.G.nodes[node]['weight'] = 1

        # Now get the node and edge positions
        node_x, node_y = self._get_node_positions()

        # Get node depths and meta
        node_name, node_hovertext, node_color = self._get_node_data()

        # Begin plotting
        edge_x, edge_y = self._get_edge_positions(self._non_skipped_edge)
        edge_trace_dark = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=self.edge_weight, color=self.dark_emphasis_color),
            hoverinfo='none',
            mode='lines',
            opacity=self.arrow_opacity,
        )
        edge_x, edge_y = self._get_edge_positions(
            lambda x: not self._non_skipped_edge(x)
        )
        edge_trace_grey = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=self.edge_weight, color=self.grey_color),
            hoverinfo='none',
            mode='lines',
            opacity=self.arrow_opacity,
        )

        layer_min, layer_max = self._get_layer_range()
        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            textposition=self.node_text_position,
            textfont_size=self.node_text_fontsize,
            hoverinfo='text',
            marker=dict(
                showscale=True,
                colorscale=self.colorscale,
                reversescale=False,
                color=[],
                size=self.node_size,
                line=dict(
                    color=self._get_border_colors(), width=self.node_border_weight
                ),
                colorbar=dict(
                    thickness=15,
                    title='Workflow Depth',
                    xanchor='left',
                    titleside='right',
                    tickmode='array',
                    tickvals=list(range(-layer_min, -layer_max - 1, -1)),
                    ticktext=list(range(layer_min, layer_max + 1, 1)),
                ),
            ),
        )
        node_trace.text = node_name
        node_trace.hovertext = node_hovertext
        node_trace.marker.color = node_color
        node_trace.marker.colorbar

        return [edge_trace_dark, edge_trace_grey, node_trace]

    def create_figure(self):
        traces = self.create_traces()
        fig = go.Figure(
            data=traces,
            layout=self._get_layout(),
        )
        fig.layout['annotations'] = self.get_annotations()
        return fig

    def _get_node_positions(self):
        positions = map(lambda n: self.G.nodes[n]['pos'], self.G.nodes)
        node_x, node_y = list(map(list, zip(*positions)))
        return node_x, node_y

    def _get_edge_positions(self, filter_fun: Callable):
        # Add position info to the edges
        edge_x = []
        edge_y = []
        for edge in list(filter(filter_fun, self.G.edges())):
            x0, y0 = self.G.nodes[edge[0]]['pos']
            x1, y1 = self.G.nodes[edge[1]]['pos']
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        return edge_x, edge_y

    def _get_node_data(self):
        node_color = list(
            map(
                lambda n: self._get_node_color(n[0], -n[1][self.node_color_key])[0],
                self.G.nodes.items(),
            )
        )
        node_name = [n for n in self.G.nodes]

        def special_labels(n):
            _, label = self._get_node_color(n, None)
            return label + '<br>' if label else ''

        node_hovertext = list(
            map(
                lambda n: (
                    special_labels(n[0])
                    + 'Stage: '
                    + str(n[0])
                    + '<br>'
                    + 'Stage Order: '
                    + str(n[1]['order'])
                    + '<br>'
                    + 'Layer: '
                    + str(n[1]['layer'])
                    + '<br>'
                    + 'Layer Order: '
                    + str(n[1]['layer_order'])
                ),
                self.G.nodes.items(),
            )
        )
        return node_name, node_hovertext, node_color

    def _get_node_color(self, n: str, default: str | int):
        if self.G.nodes[n]['skip_stages']:
            return '#5A5A5A', 'Skip stage'
        elif self.G.nodes[n]['only_stages']:
            return '#5053f8', 'Only run this stage'
        elif self.G.nodes[n]['first_stages']:
            return '#33c584', 'First stage'
        elif self.G.nodes[n]['last_stages']:
            return '#e93e2e', 'Last stage'
        elif self.G.nodes[n]['skipped']:
            return self.grey_color, 'Skipped'
        else:
            return default, None

    def _get_border_color(self, n: str):
        if self.G.nodes[n]['skipped']:
            return self.grey_color
        else:
            return self.dark_emphasis_color

    def _non_skipped_edge(self, edge) -> bool:
        return (
            not self.G.nodes[edge[0]]['skipped']
            and not self.G.nodes[edge[1]]['skipped']
        )

    def _non_skipped_node(self, node: str) -> bool:
        return not self.G.nodes[node]['skipped']

    def _get_edge_color(self, edge):
        return (
            self.dark_emphasis_color
            if self._non_skipped_edge(edge)
            else self.grey_color
        )

    def _get_node_colors(self):
        return [self._get_node_color(n)[0] for n in self.G.nodes]

    def _get_border_colors(self):
        return [self._get_border_color(n) for n in self.G.nodes]

    def _get_edge_colors(self):
        return [self._get_edge_color(e) for e in self.G.edges]

    def _recalculate_depth(self, new_key: str):
        for layer, nodes in enumerate(nx.topological_generations(self.G)):
            for node in nodes:
                self.G.nodes[node][new_key] = layer

    def _calculate_depth_order(self, layer_key: str, new_key: str):
        # Add position info to the graph nodes
        pos = nx.multipartite_layout(
            self.G,
            subset_key=layer_key,
            align=self.align,
            scale=self.layout_scale,
        )
        for n, p in pos.items():
            self.G.nodes[n]['pos'] = p

        # Get all the node meta, add the name as well so we can just pass around values
        nodes = dict(self.G.nodes.items())
        nodes = {n: dict(meta, name=n) for n, meta in nodes.items()}

        # Group by partite_key
        sorted_nodes = sorted(nodes.values(), key=lambda n: n[layer_key])
        by_depth = groupby(sorted_nodes, lambda x: x[layer_key])

        # Go through each layer group
        for _, group in by_depth:
            # Extract nodes and node positions (sorted)
            group_nodes = list(group)
            group_pos = sorted([n['pos'] for n in group_nodes], key=lambda p: p[0])

            # Do a layer sort
            nodes = self._node_layer_sort(group_nodes)

            # Iterate through all running jobs, then skipped
            # Set layer_order to it's order index in that layer
            for depth_order, node in enumerate(nodes):
                idx = depth_order
                self.G.nodes[node['name']][new_key] = idx
                self.G.nodes[node['name']]['pos'] = group_pos[idx]

    def get_annotations(self, xref='x', yref='y'):
        def pts(edge, key):
            pts = {
                'p1.x': self.G.nodes[edge[0]]['pos'][0],
                'p1.y': self.G.nodes[edge[0]]['pos'][1],
                'p2.x': self.G.nodes[edge[1]]['pos'][0],
                'p2.y': self.G.nodes[edge[1]]['pos'][1],
            }
            return pts.get(key)

        return [
            go.layout.Annotation(
                ax=(pts(edge, 'p1.x') * 7 + pts(edge, 'p2.x') * 1) / 8,
                ay=(pts(edge, 'p1.y') * 7 + pts(edge, 'p2.y') * 1) / 8,
                x=(pts(edge, 'p1.x') * 5 + pts(edge, 'p2.x') * 3) / 8,
                y=(pts(edge, 'p1.y') * 5 + pts(edge, 'p2.y') * 3) / 8,
                axref=xref,
                ayref=yref,
                xref=xref,
                yref=yref,
                showarrow=True,
                arrowhead=3,
                arrowsize=4,
                arrowwidth=self.edge_weight,
                arrowcolor=self._get_edge_color(edge),
                opacity=self.arrow_opacity,
            )
            for edge in self.G.edges
        ]

    def _node_layer_sort(self, nodes):
        layer_num_parity = nodes[0][self.partite_key] % 2
        degree = self.G.out_degree([n['name'] for n in nodes])
        deg_order = sorted(
            nodes, key=lambda n: degree[n['name']], reverse=layer_num_parity
        )

        # Set largest degree towards the middle
        return deg_order[len(deg_order) % 2 :: 2] + deg_order[::-2]

        # Sort all non-skipped stages
        sorted_group = list(
            sorted(
                filter(lambda x: not x['skipped'], nodes),
                key=lambda x: x['order'],
            )
        )

        # Then sort all skipped nodes
        skipped_group = list(
            sorted(
                filter(lambda x: x['skipped'], nodes),
                key=lambda x: x['order'],
            )
        )
        return sorted_group + skipped_group

    def _get_layer_range(self):
        nodes = set(filter(self._non_skipped_node, self.G.nodes))
        layer = nx.get_node_attributes(self.G, self.partite_key)
        layer = set([val for node, val in layer.items() if node in nodes])
        return min(layer), max(layer)

    def _get_node_file_link(self, node: str):
        def normalize(string: str):
            return string.lower().replace('-', '').replace('_', '').replace('.py', '')

        files = os.scandir(self.stages_folder)
        files = {normalize(os.path.basename(f)): f for f in files}
        print(files)
        