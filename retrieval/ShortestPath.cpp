// lemon
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/maps.h>
#include <lemon/dijkstra.h>

using namespace std;
using namespace lemon;

int main()
{
    ListDigraph g;
    ListDigraph::NodeMap<string> cuis(g);
    ListDigraph::ArcMap<float> arcWeights(g);
    
    ListDigraph::Node s = g.addNode();
    ListDigraph::Node t = g.addNode();
    ListDigraph::Node x = g.addNode();
    
    ListDigraph::Arc a = g.addArc(s,t);
    ListDigraph::Arc a2 = g.addArc(t,x);
    
    arcWeights[a] = 2;
    arcWeights[a2] = 1;
    
    float d = 0.5;
    
    bool reached = dijkstra(g,arcWeights).dist(d).run(s,x);
    
    cout << d << endl;
    
    return 0;
}