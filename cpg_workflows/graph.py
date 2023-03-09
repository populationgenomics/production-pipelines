"""
Code to display the prod pipelines graph

Heavily inspired by:
https://towardsdatascience.com/visualize-hierarchical-data-using-plotly-and-datapane-7e5abe2686e1
"""

import networkx as nx
import plotly.graph_objs as go

DEFAULT_GML_OUTPUT = './.workflow.gml'


class GraphPlot:
    def __init__(self, G: nx.DiGraph):
        self.G = G
        self.title = '8-Level network graph'
        self.title_fontsize = 16
        self.node_size = 10
        self.node_border_weight = 2
        self.edge_weight = 0.5
        self.edge_color = '#888'

    def display_graph(self):
        # Add position info to the graph nodes
        pos = nx.kamada_kawai_layout(self.G, weight=None)
        for n, p in pos.items():
            self.G.nodes[n]['pos'] = p

        # Add weight and depth attributes to the nodes
        for node in self.G.nodes:
            self.G.nodes[node]['weight'] = 1
            self.G.nodes[node]['depth'] = 1

        # Now get the node and edge positions
        node_x, node_y = self._get_node_positions()
        edge_x, edge_y = self._get_edge_positions()

        # Get node depths and meta
        node_depth, node_text = self._get_node_data()

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
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=self.node_size,
                line_width=self.node_border_weight,
                colorbar=dict(
                    thickness=15,
                    title='Indent Level',
                    xanchor='left',
                    titleside='right',
                ),
            ),
        )
        node_trace.marker.color = node_depth
        node_trace.text = node_text

        # Draw the graph
        graph = go.Figure(
            data=[edge_trace, node_trace],
            layout=go.Layout(
                title=self.title,
                titlefont_size=self.title_fontsize,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
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
        node_depth = list(map(lambda n: n[1]['depth'], self.G.nodes.items()))
        node_text = list(
            map(
                lambda n: '{0} {1:.2f}%'.format(n[0], n[1]['weight']),
                self.G.nodes.items(),
            )
        )
        return node_depth, node_text
