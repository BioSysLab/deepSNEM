#!/usr/bin/env python
# -*- coding: utf-8 -*-

from models.deep_graph_infomax.infomax import SNInfomax
from models.graph_transformer.euclidean_graph_transformer import GraphTransformerEncoder

encoder = GraphTransformerEncoder(1, 4, 512, None, None, True)
sninfomax = SNInfomax(512, encoder, summary=None)

if __name__ == "__main__":
    print(list(sninfomax.modules()))
