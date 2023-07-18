class Louvain {
  static int m = 28;
  List<Node> nodes;
  List<Community> communities;

  void exec(int times) {
    for (int i = 0; i < times; i++) {
      modularityOptimization();
      communityAggregation();
    }
  }

  void modularityOptimization() {
    bool move;
    do {
      move = false;
      for (Node node in nodes) {
        Node? tempNode;
        double max_delta = 0;

        for (Node adj in node.next) {
          if (node.community.id != adj.community.id) {
            double takeOut = node.community.modularityOut(node);
            double putIn = adj.community.modularityIn(node);
            double delta = takeOut + putIn;

            if (delta > max_delta) {
              tempNode = adj;
              max_delta = delta;
            }
          }
        }

        if (tempNode != null) {
          move = true;
          tempNode.community.nodes.add(node);
          node.community.nodes.remove(node);
          node.community = tempNode.community;
        }
      }
    } while (move);
  }

  void communityAggregation() {
    Map<int, int> index = <int, int>{};
    List<Node> new_nodes = [];
    List<Community> new_communities = [];
    int id = 0;

    for (Node node in nodes) {
      if (!index.containsKey(node.community.id)) {
        index[node.community.id] = id;
        new_nodes.add(Node(id, node.community.sumIn() + node.self_weight));
        id++;
      } else {
        new_nodes[index[node.community.id]!].self_weight += node.self_weight;
      }
    }

    for (Node node in nodes) {
      for (Node adj in node.next) {
        if (node.community.id != adj.community.id) {
          Node curr = new_nodes[index[node.community.id]!];
          Node next = new_nodes[index[adj.community.id]!];

          if (!curr.weight.containsKey(next)) {
            curr.add(next, node.weight[adj]!);
          } else {
            curr.weight[next] = curr.weight[next]! + node.weight[adj]!;
          }
        }
      }
    }

    for (Node node in new_nodes) {
      new_communities.add(node.community);
    }

    nodes = new_nodes;
    communities = new_communities;
  }

  Louvain(this.nodes, this.communities);
}

class Node {
  int id;
  int self_weight;
  List<Node> next = [];
  Map<Node, int> weight = <Node, int>{};
  late Community community;

  Node(this.id, this.self_weight) {
    community = Community(this);
  }

  void add(Node node, int value) {
    next.add(node);
    weight[node] = value;
  }

  int kin(Community community) {
    int result = 0;

    for (Node node in community.nodes) {
      if (next.contains(node)) {
        result += weight[node]!;
      }
    }

    return result * 2;
  }

  int ki() {
    int result = 0;

    for (Node node in next) {
      result += weight[node]!;
    }

    return result;
  }

  double modularity() {
    return (self_weight / 2 * Louvain.m) - (ki() / (2 * Louvain.m)) * (ki() / (2 * Louvain.m));
  }

  int getSumWeight() {
    int result = 0;
    for (Node node in next) {
      result += weight[node]!;
    }

    return result;
  }
}

class Community {
  late int id;
  List<Node> nodes = [];

  Community(Node node) {
    id = node.id;
    nodes.add(node);
  }

  void add(Node node) {
    node.community = this;

    nodes.add(node);
  }

  int sumIn() {
    int result = 0;

    for (Node node in nodes) {
      for (Node adj in node.next) {
        if (nodes.contains(adj)) {
          result += node.weight[adj]!;
        }
      }
    }

    return result;
  }

  int sumTot() {
    int result = 0;

    for (Node node in nodes) {
      result += node.getSumWeight();
    }

    return result;
  }

  double modularity() {
    return (sumIn() / (2 * Louvain.m)) -
        (sumTot() / (2 * Louvain.m)) * (sumTot() / (2 * Louvain.m));
  }

  double modularityIn(Node node) {
    double before = modularity() + node.modularity();
    double after = ((sumIn() + node.kin(this)) / (2 * Louvain.m)) -
        ((sumTot() + node.ki()) / (2 * Louvain.m)) * ((sumTot() + node.ki()) / (2 * Louvain.m));
    return after - before;
  }

  double modularityOut(Node node) {
    nodes.remove(node);
    double after = modularity() + node.modularity();
    double before = ((sumIn() + node.kin(this)) / (2 * Louvain.m)) -
        ((sumTot() + node.ki()) / (2 * Louvain.m)) * ((sumTot() + node.ki()) / (2 * Louvain.m));
    nodes.add(node);
    return after - before;
  }
}
