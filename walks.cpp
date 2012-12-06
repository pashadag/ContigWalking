#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;
ostringstream oMsg;
string sbuf;
#include "../Tools/dbtools/defs.h"
#include "../Tools/dbtools/graph.h"
#include "../Tools/dbtools/collapsed_db_from_genome.h"

typedef set<  int > adjList; // reference into edges
typedef vector< adjList > adjType;

int numNodes;

spAdjType convertGraphForSp(const collapsed_db_from_genome & db, const adjType & g) {
	spAdjType h(g.size());
	for (int i = 0; i < g.size(); i++) {
		for (adjList::const_iterator j = g[i].begin(); j != g[i].end(); j++) {
			// this could be an out-neighbor graph or an in-neighbor graph;
			int connectedNode = db.edges[*j].to;
			if (connectedNode == i) connectedNode = db.edges[*j].from;
			h[i].insert(make_pair(connectedNode, db.edges[*j].length));
		}
	}
	return h;
}
void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " [parameters]" << endl;
	cerr << "Program descrption.\n";
	exit(1);
}

void renumberEdges(collapsed_db_from_genome & db) {  //remove isolated vertices
	cout << "Before renumbering: " << numNodes << " nodes, " << db.edges.size() << " edges.\n";

	vector<int> old2new(numNodes, -1);
	for (int i = 0; i < db.edges.size(); i++) {
		int from = db.edges[i].from;
		int to   = db.edges[i].to;
		old2new[from] = -2;
		old2new[to]   = -2;
	}

	int curNum = 0;
	for (int i = 0; i < numNodes; i++) {
		if (old2new[i] != -1) {
			old2new[i] = curNum++;
		}
	}
	
	for (int i = 0; i < db.edges.size(); i++) {
		db.edges[i].from = old2new.at(db.edges[i].from);
		db.edges[i].to   = old2new.at(db.edges[i].to);
	}

	numNodes = curNum;

	cout << "After renumbering: " << numNodes << " nodes, " << db.edges.size() << " edges.\n";
}


adjType getOutGraph(collapsed_db_from_genome & db, int node = -1) {  //optionally removes node from the graph
	adjType graph;
	graph.resize(numNodes);
	for (int i = 0; i < db.edges.size(); i++) {
		if (db.edges[i].from != node && db.edges[i].to != node) {
			graph[db.edges[i].from].insert(i);
		}
	}
	return graph;
}

adjType getInGraph(collapsed_db_from_genome & db, int node = -1) {
	adjType graph;
	graph.resize(numNodes);
	for (int i = 0; i < db.edges.size(); i++) {
		if (db.edges[i].from != node && db.edges[i].to != node) {
			graph[db.edges[i].to].insert(i);
		}
	}
	return graph;
}


void debugTest( collapsed_db_from_genome & db) {
	for (int i = 0; i < db.edges.size(); i++) {
		edge e = db.edges[i];
		if (e.from == 23 && e.to == 27) cout << "23-->27 " << db.genome.substr(e.seqStart, e.length) << endl;
		if (e.from == 27 && e.to == 10) cout << "27-->10 " << db.genome.substr(e.seqStart, e.length) << endl;
		if (e.from == 27 && e.to == 11) cout << "27-->11 " << db.genome.substr(e.seqStart, e.length) << endl;
		if (e.from == 11 && e.to == 10) cout << "11-->10 " << db.genome.substr(e.seqStart, e.length) << endl;
	}
}

set<pair<int,int> > dupEdges;

bool isEdgeMultiple(collapsed_db_from_genome & db, pair<int,int> target) {
	if (dupEdges.size() == 0) {
		//first time around, lets build dupEdges datastructure 
		for (int i = 1; i < db.edges.size(); i++) {
			if (db.edges[i].from == db.edges[i-1].from && db.edges[i].to  == db.edges[i-1].to ) {
				dupEdges.insert(make_pair(db.edges[i].from, db.edges[i].to));
			}
		}
	}
	set<pair<int,int> >::iterator it = dupEdges.find(target);
	return (it != dupEdges.end());
}

void process_node(collapsed_db_from_genome & db, adjType & outGraph, adjType & inGraph, long mainNode) {

	//Recall:
	//    typedef set< int > adjList;
	//    typedef vector< adjList > adjType;

	//if there is a loop at this node, abandon
	for (adjList::iterator idxOutEdge = outGraph[mainNode].begin(); idxOutEdge != outGraph[mainNode].end(); idxOutEdge++) {
		if (db.edges[*idxOutEdge].to == mainNode) {
			cout << "Node " << mainNode << " has loop, abandoning.\n";
			return;
		}
	}

	
	set<pair<int,int> > reachable;
	//vector<bool> dummy(numNodes, false);
	//vector< vector<bool> > reachable (numNodes, dummy);

	//remove mainNode from the graph
	spAdjType spModOutGraph = convertGraphForSp(db, getOutGraph(db, mainNode));

	//fill in reachability matrix
	for (adjList::iterator idxOutEdge = outGraph[mainNode].begin(); idxOutEdge != outGraph[mainNode].end(); idxOutEdge++) {
		int outNode = db.edges[*idxOutEdge].to;
		set<int> found;
		vector<int> dist;
		sp(outNode, dist, spModOutGraph, found, -1); 
		for (adjList::iterator idxInEdge = inGraph[mainNode].begin(); idxInEdge != inGraph[mainNode].end(); idxInEdge++) {
			int inNode = db.edges[*idxInEdge].from;
			if (dist[inNode] != -1) {
				reachable.insert(make_pair(outNode,inNode));
				//reachable[outNode][inNode] = true;
			}
		}
	}

	//cout << print_sparce(reachable, true) << endl;
	

	//connect edges according to reachability
	for (adjList::iterator idxOutEdge = outGraph[mainNode].begin(); idxOutEdge != outGraph[mainNode].end(); idxOutEdge++) {
		int outNode = db.edges[*idxOutEdge].to;
		if (isEdgeMultiple(db, make_pair(mainNode, outNode))) continue;
		for (adjList::iterator idxInEdge = inGraph[mainNode].begin(); idxInEdge != inGraph[mainNode].end(); idxInEdge++) {
			int inNode = db.edges[*idxInEdge].from;
			bool contig = true;
			if (isEdgeMultiple(db, make_pair(inNode, mainNode))) continue;

			for (adjList::iterator idxOutEdge2 = outGraph[mainNode].begin(); idxOutEdge2 != outGraph[mainNode].end(); idxOutEdge2++) {
				int outNode2 = db.edges[*idxOutEdge2].to;
				for (adjList::iterator idxInEdge2 = inGraph[mainNode].begin(); idxInEdge2 != inGraph[mainNode].end(); idxInEdge2++) {
					int inNode2 = db.edges[*idxInEdge2].from;

					if (inNode != inNode2 && outNode != outNode2) {
						set<pair<int,int> >::iterator it = reachable.find(make_pair(outNode2, inNode2));
						if (it != reachable.end()) {
						//if (reachable[outNode2][inNode2]) {
							contig = false;
							goto NEXT;
						}
					}
				}
			}
NEXT:   	
			if (contig) {
				cout << "Join: " << inNode << " -> " << mainNode << " -> " << outNode << " with lengths " << db.edges[*idxInEdge].length << " and " << db.edges[*idxOutEdge].length << endl;
				edge * inEdge  = &db.edges[*idxInEdge];
				edge * outEdge = &db.edges[*idxOutEdge];
				cout << "CONTIG " << db.genome.substr(inEdge->seqStart, inEdge->length) << db.genome.substr(outEdge->seqStart, outEdge->length) << endl;
			}
		}
	}

}


void removeSplitJoinNodes(collapsed_db_from_genome & db) {

	cout << "Before split/join removal: " << numNodes << " nodes, " << db.edges.size() << " edges.\n";
	adjType outGraph = getOutGraph(db);
	adjType inGraph  = getInGraph(db);
	int nodesRemoved = 0;

	//remove split nodes
	for (int i = 0; i < numNodes; i++) {
		if (inGraph[i].size() != 1)  continue;
		int idxInEdge = *inGraph[i].begin();
		edge * inEdge = &db.edges[idxInEdge];

		if (outGraph[i].size() == 1) continue; //this might happen if its a vertex with just a loop

		for (adjList::iterator j = outGraph[i].begin(); j != outGraph[i].end(); j++) {
			//make a new edge
			edge * outEdge = &db.edges[*j];
			outEdge->from = inEdge->from;
			outEdge->length += inEdge->length;
			outEdge->seqStart -= inEdge->length;

			outGraph[inEdge->from].insert(*j);
		}
		outGraph[inEdge->from].erase(outGraph[inEdge->from].find(idxInEdge));
		outGraph[i].clear();
		inGraph[i].clear();
		inEdge->from = -1; //mark inEdge for removal

		nodesRemoved++;
	}

	//remove join nodes
	//this is cut and paste with symmetrical changes (from/to, in/out).  The only other diff is that seqStart no longer needs to be changed
	for (int i = 0; i < numNodes; i++) {
		if (outGraph[i].size() != 1)  continue;
		int idxOutEdge = *outGraph[i].begin();
		edge * outEdge = &db.edges[idxOutEdge];

		if (inGraph[i].size() == 1) continue; //this might happen if its a vertex with just a loop

		for (adjList::iterator j = inGraph[i].begin(); j != inGraph[i].end(); j++) {
			//make a new edge
			edge * inEdge = &db.edges[*j];
			inEdge->to = outEdge->to;
			inEdge->length += outEdge->length;
			//inEdge->seqStart -= inEdge->length;

			inGraph[outEdge->to].insert(*j);
		}
		inGraph[outEdge->to].erase(inGraph[outEdge->to].find(idxOutEdge));
		inGraph[i].clear();
		outGraph[i].clear();
		outEdge->to = -1; //mark inEdge for removal

		nodesRemoved++;
	}


	//rebuild edges
	vector<edge> newEdges;
	for (int i = 0; i < db.edges.size(); i++) {
		if (db.edges[i].to != -1 && db.edges[i].from != -1) newEdges.push_back(db.edges[i]);
	}
	db.edges = newEdges;

	cout << "After split/join removal: " << numNodes - nodesRemoved << " nodes, " << db.edges.size() << " edges.\n";
}


int main(int argc, char * argv[]) {


	assert(argc == 3);
	string genomeBase = argv[1];
	int kmersize = atoi(argv[2]);
	cout << "Reading genome...\n";
	string genome = read_genome(genomeBase + ".fa");

	cout << "Loading repeat graph..\n";
	collapsed_db_from_genome db(genome, kmersize);
	//db.build();
	db.load(genomeBase + ".db");
	for (int i = 0; i < db.edges.size(); i++) numNodes = max(max(numNodes, db.edges[i].from + 1), db.edges[i].to + 1);

	cout << "Removing identical edges in repeat graph..\n";
	db.removeIdenticalEdges();
	//debugTest(db);

	cout << "Writing out dot graph...\n";
	db.writeDOT(genomeBase + ".iden.dot");

	cout << "Writing out edge sequences...\n";
	db.writeEdgeSeqs(genomeBase + ".iden.contigs");

	cout << "Removing split nodes...\n";
	removeSplitJoinNodes(db);
	renumberEdges(db);
	sort(db.edges.begin(), db.edges.end());

	cout << "Writing out edge sequences...\n";
	db.writeEdgeSeqs(genomeBase + ".split.contigs");

	cout << "Writing out dot graph...\n";
	db.writeDOT(genomeBase + ".split.dot");


	adjType inGraph, outGraph;
	inGraph = getInGraph(db);
	outGraph = getOutGraph(db);
	for (int i = 0; i < numNodes; i ++) {
		cout << "Processing node " << i << "..." << endl;
		process_node(db, outGraph, inGraph, i);
	}




	return 0;
}


