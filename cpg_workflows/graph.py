import networkx as nx
import plotly.graph_objs as go

DEFAULT_GML_OUTPUT = './.workflow.gml'


class GraphPlot:
    def __init__(self, G: nx.DiGraph):
        self.G = G
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

        # Now get the node and edge positions
        node_x, node_y = self._get_node_positions()
        edge_x, edge_y = self._get_edge_positions()

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

    def _get_node_positions(self):
        positions = map(lambda n: n['pos'], self.G.nodes)
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
