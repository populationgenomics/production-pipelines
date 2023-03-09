"""
Code to display the prod pipelines graph

Heavily inspired by:
https://towardsdatascience.com/visualize-hierarchical-data-using-plotly-and-datapane-7e5abe2686e1
"""

from itertools import groupby

import networkx as nx
import plotly.graph_objs as go

DEFAULT_GML_OUTPUT = './.workflow.gml'


class GraphPlot:
    def __init__(self, G: nx.DiGraph):
        self.G = G

        # Text and sizes
        self.title = 'Workflow Graph'
        self.title_fontsize = 16
        self.node_size = 50
        self.node_border_weight = 2
        self.edge_weight = 1

        # Colors
        self.colorscale = 'Blugrn'
        self.skipped_color = '#D3D3D3'
        self.edge_color = '#888'
        self.default_border_color = '#888'

        # Layout
        self.align = 'horizontal'

        # Graph transforms

        # Add name as attribute
        nx.set_node_attributes(self.G, {n: n for n in self.G.nodes}, name='name')

        # We invert the depths so that the starting stages have depth 0
        self._invert_depth()

        # Calculate the depth_order
        self._calculate_depth_order()

        # Create a new attribute for the position:
        # depth.order
        # self.partite_key = 'depth'
        # nx.set_node_attributes(
        #     self.G,
        #     {
        #         n[0]: float(f'{n[1]["depth"]}.{n[1]["order"]}')
        #         for n in self.G.nodes.items()
        #     },
        #     name=self.partite_key,
        # )

    def display_graph(self):
        # Add position info to the graph nodes

        for n, meta in self.G.nodes.items():
            self.G.nodes[n]['pos'] = (meta['depth_order'], meta['depth'])

        # pos = nx.multipartite_layout(
        #     self.G, subset_key=self.partite_key, align=self.align
        # )
        # pos = nx.kamada_kawai_layout(
        #     self.G,
        # )
        # for n, p in pos.items():
        #     self.G.nodes[n]['pos'] = p

        # Add weight and depth attributes to the nodes
        for node in self.G.nodes:
            self.G.nodes[node]['weight'] = 1

        # Now get the node and edge positions
        node_x, node_y = self._get_node_positions()
        edge_x, edge_y = self._get_edge_positions()

        # Get node depths and meta
        node_name, node_hovertext, node_color = self._get_node_data()

        # Begin plotting
        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=self.edge_weight, color=self.edge_color),
            hoverinfo='none',
            mode='lines',
        )

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            textposition="bottom center",
            hoverinfo='text',
            marker=dict(
                showscale=True,
                colorscale=self.colorscale,
                reversescale=True,
                color=[],
                size=self.node_size,
                line_width=self.node_border_weight,
                colorbar=dict(
                    thickness=15,
                    title='Workflow Depth',
                    xanchor='left',
                    titleside='right',
                ),
            ),
        )
        node_trace.text = node_name
        node_trace.hovertext = node_hovertext
        node_trace.marker.color = node_color

        # Draw the graph
        annotations = self._get_annotations()

        graph = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title=self.title,
                titlefont_size=self.title_fontsize,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(
                    showgrid=False,
                    zeroline=False,
                    showticklabels=False,
                    autorange='reversed',
                ),
                annotations=annotations,
            ),
        )

        graph.show()

    def _get_node_positions(self):
        positions = map(lambda n: self.G.nodes[n]['pos'], self.G.nodes)
        node_x, node_y = list(map(list, zip(*positions)))
        return node_x, node_y

    def _get_edge_positions(self):
        # Add position info to the edges
        edge_x = []
        edge_y = []
        for edge in self.G.edges():
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
                lambda n: self.skipped_color if n[1]['skipped'] else n[1]['depth'],
                self.G.nodes.items(),
            )
        )
        node_name = [n for n in self.G.nodes]
        node_hovertext = list(
            map(
                lambda n: (
                    'Stage: '
                    + str(n[0])
                    + '<br>'
                    + 'Stage Order: '
                    + str(n[1]['order'])
                    + '<br>'
                    + 'Depth: '
                    + str(n[1]['depth'])
                    + '<br>'
                    + 'Depth Order: '
                    + str(n[1]['depth_order'])
                ),
                self.G.nodes.items(),
            )
        )
        return node_name, node_hovertext, node_color

    def _get_border_colors(self):
        def get_border_color(n):
            if self.G.nodes[n]['skip_stages']:
                return '#5A5A5A'
            elif self.G.nodes[n]['only_stages']:
                return '#00008B'
            elif self.G.nodes[n]['first_stages']:
                return '#023020'
            elif self.G.nodes[n]['last_stages']:
                return '#8B0000'
            else:
                return self.default_border_color

        return {n: get_border_color(n) for n in self.G.nodes}

    def _invert_depth(self):
        depths = nx.get_node_attributes(self.G, 'depth')
        max_depth = max(depths.values())
        new_depths = {stage: abs(max_depth - depth) for stage, depth in depths.items()}
        nx.set_node_attributes(self.G, new_depths, name='depth')

    def _calculate_depth_order(self):
        nodes = dict(self.G.nodes.items())
        nodes = {n: dict(meta, name=n) for n, meta in nodes.items()}

        by_depth = groupby(nodes.values(), lambda x: x['depth'])

        for _, group in by_depth:
            for depth_order, node in enumerate(sorted(group, key=lambda x: x['order'])):
                self.G.nodes[node['name']]['depth_order'] = depth_order

    def _get_annotations(self):
        return [
            dict(
                ax=(self.G.nodes[edge[0]]['pos'][0] + self.G.nodes[edge[1]]['pos'][0])
                / 2,
                ay=(self.G.nodes[edge[0]]['pos'][1] + self.G.nodes[edge[1]]['pos'][1])
                / 2,
                axref='x',
                ayref='y',
                x=(
                    self.G.nodes[edge[1]]['pos'][0] * 3
                    + self.G.nodes[edge[0]]['pos'][0]
                )
                / 4,
                y=(
                    self.G.nodes[edge[1]]['pos'][1] * 3
                    + self.G.nodes[edge[0]]['pos'][1]
                )
                / 4,
                xref='x',
                yref='y',
                showarrow=True,
                arrowhead=3,
                arrowsize=4,
                arrowwidth=self.edge_weight,
                arrowcolor=self.edge_color,
                opacity=1,
            )
            for edge in self.G.edges
        ]
