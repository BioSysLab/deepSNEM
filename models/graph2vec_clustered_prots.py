"""Graph2Vec module."""

import json
import glob
import hashlib
import pandas as pd
import networkx as nx
from tqdm import tqdm
from joblib import Parallel, delayed
from param_parser import parameter_parser
from gensim.models.doc2vec import Doc2Vec, TaggedDocument

class WeisfeilerLehmanMachine:
    """
    Weisfeiler Lehman feature extractor class.
    """
    def __init__(self, graph, features, iterations):
        """
        Initialization method which also executes feature extraction.
        :param graph: The Nx graph object.
        :param features: Feature hash table.
        :param iterations: Number of WL iterations.
        """
        self.iterations = iterations
        self.graph = graph
        self.features = features
        self.nodes = self.graph.nodes()
        self.extracted_features = [str(v) for k, v in features.items()]
        self.do_recursions()

    def do_a_recursion(self):
        """
        The method does a single WL recursion.
        :return new_features: The hash table with extracted WL features.
        """
        new_features = {}
        for node in self.nodes:
            nebs = self.graph.neighbors(node)
            degs = [self.features[neb] for neb in nebs]
            features = [str(self.features[node])]+sorted([str(deg) for deg in degs])
            features = "_".join(features)
            hash_object = hashlib.md5(features.encode())
            hashing = hash_object.hexdigest()
            new_features[node] = hashing
        self.extracted_features = self.extracted_features + list(new_features.values())
        return new_features

    def do_recursions(self):
        """
        The method does a series of WL recursions.
        """
        for _ in range(self.iterations):
            self.features = self.do_a_recursion()

def net2json(df):
    nodes=list(set(list(df['node1']) + list(df['node2'])))
    #read pickle in the folder you are
    #import pickle
    #filep=open('ReSimNet-Dataset.pkl','rb')
    #prot_enc=pickle.load(filep)
    #filep.close()
    prot_enc=pd.read_csv("prot_clustering.csv",index_col=0)
    new=[]
    cl_new=[]
    cl=1500
    if "Perturbation" in nodes:
        nodes.remove("Perturbation")
    #nodes[nodes.index("Perturbation")]=nodes[0]
    #nodes[0]="Perturbation"
    feats={}
    for j,x in enumerate(nodes):
        if (x in list(df['node1'])):
            activity=list(df['activity1'])[list(df['node1']).index(x)]
        else:
            activity=list(df['activity2'])[list(df['node2']).index(x)]
        sign = lambda a: '-' if (int(a)==-1) else '+'
        if (x in list(prot_enc['protein'])):
            pr=list(prot_enc['cluster'])[list(prot_enc['protein']).index(x)]
        else:
            if (x in new):
                pr=cl_new[new.index(x)]
            else:
                new.append(x)
                cl+=1
                cl_new='Cluster'+str(cl)
                pr=cl_new[new.index(x)]
        feats.update({"%s"%j:pr+sign(activity)})
    net={}
    net.update({"edges":[]})
    net.update({"features":feats})
    for j in range(1,len(df)):
        if ((df['node1'][j]!="Perturbation") and (df['node2'][j]!="Perturbation")):
            ind1=nodes.index(df['node1'][j])
            ind2=nodes.index(df['node2'][j])
            net["edges"].append([ind1,ind2])
    return(net)

def dataset_reader(path):
    """
    Function to read the graph and features from a json file.
    :param path: The path to the graph json.
    :return graph: The graph object.
    :return features: Features hash table.
    :return name: Name of the graph.
    """
    emb=path.strip(".csv").split("/")[-1]
    emb=emb.split("_")[-1]
    name = path.strip(".csv").split("/")[-2]+"_emb_"+str(emb)
    data = pd.read_csv(path,index_col=0)
    data=net2json(data)
    graph = nx.from_edgelist(data["edges"])

    if "features" in data.keys():
        features = data["features"]
    else:
        features = nx.degree(graph)

    features = {int(k): v for k, v in features.items()}
    return graph, features, name

def feature_extractor(path, rounds):
    """
    Function to extract WL features from a graph.
    :param path: The path to the graph json.
    :param rounds: Number of WL iterations.
    :return doc: Document collection object.
    """
    graph, features, name = dataset_reader(path)
    machine = WeisfeilerLehmanMachine(graph, features, rounds)
    doc = TaggedDocument(words=machine.extracted_features, tags=["g_" + name])
    return doc

def save_embedding(output_path, model, files, dimensions):
    """
    Function to save the embedding.
    :param output_path: Path to the embedding csv.
    :param model: The embedding model object.
    :param files: The list of files.
    :param dimensions: The embedding dimension parameter.
    """
    out = []
    for f in files:
        #identifier = f.split("/")[-1].strip(".json")
        emb=f.strip(".csv").split("/")[-1]
        emb=emb.split("_")[-1]
        identifier = f.split("/")[-2]+"_emb_"+str(emb)
        out.append([identifier] + list(model.docvecs["g_"+identifier]))
    column_names = ["type"]+["x_"+str(dim) for dim in range(dimensions)]
    out = pd.DataFrame(out, columns=column_names)
    out = out.sort_values(["type"])
    out.to_csv(output_path, index=None)

def main(args):
    """
    Main function to read the graph list, extract features.
    Learn the embedding and save it.
    :param args: Object with the arguments.
    """
    #graphs = glob.glob(args.input_path + "*.json")
    in_graphs = args.input_path
    graphs=pd.read_csv(in_graphs,index_col=0)
    graphs= list(graphs['path_list'])
    print("\nFeature extraction started.\n")
    #document_collections = Parallel(n_jobs=args.workers)(delayed(feature_extractor)(g, args.wl_iterations) for g in tqdm(graphs))
    document_collections = Parallel(n_jobs=args.workers)(delayed(feature_extractor)(g, args.wl_iterations) for g in tqdm(graphs))
    print("\nOptimization started.\n")

    model = Doc2Vec(document_collections,
                    vector_size=args.dimensions,
                    window=0,
                    min_count=args.min_count,
                    dm=0,
                    sample=args.down_sampling,
                    workers=args.workers,
                    epochs=args.epochs,
                    alpha=args.learning_rate)

    save_embedding(args.output_path, model, graphs, args.dimensions)

if __name__ == "__main__":
    args = parameter_parser()
    main(args)
